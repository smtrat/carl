/* 
 * File:   TermAdditionManager.h
 * Author: Florian Corzilius
 *
 * Created on October 30, 2014, 7:20 AM
 */

#pragma once 
#include "../core/Term.h"
#include <vector>
#include <unordered_map>
#include <mutex>

namespace carl
{

/**
 * Class to manage term addition.
 */
template<typename Polynomial>
class TermAdditionManager
{
    private:
        
        /// A hash map of shared pointer of monomials to shared pointer of terms. 
        typedef std::unordered_map<std::shared_ptr<const Monomial>,std::shared_ptr<const Term<typename Polynomial::CoeffType>>, std::equal_to<std::shared_ptr<const Monomial>>> MapType;
    
        /// Id of the next free map.
        size_t mNextMapId;
        /** 
         * If the ith map is used by the polynomial p, this map stores a reference to p at the ith entry of this vector. 
         * This is only for asserting that no illegal access occurs.
         */
        std::vector<const Polynomial*> mUsers;
        /// Stores the maps for the detection of equal terms.
        std::vector<MapType> mTermMaps;
        /** 
         * Stores the leading term of the currently added terms for each map. It points to the end of its corresponding
         * map, if no leading term has been found (a constant term will not be taken as a leading term).
         */
        std::vector<typename MapType::iterator> mLeadingTerms;
        /** 
         * Stores the constant term of the currently added terms for each map. It points to the end of its corresponding
         * map, if the constant term is currently 0.
         */
        std::vector<typename MapType::iterator> mConstantTerms;
        /// Mutex for allocation of the map ids.
        std::mutex mMutex;
        
    public:
            
        /**
         * Constructs a term addition manager.
         * @param _numberOfMaps The number of equal-term-detection-maps this manager manages.
         */
        TermAdditionManager( size_t _numberOfMaps = 1 ):
            mNextMapId( 0 ),
            mUsers( _numberOfMaps, nullptr ),
            mTermMaps( _numberOfMaps, MapType() ),
            mLeadingTerms( _numberOfMaps, mTermMaps.end() ),
            mConstantTerms( _numberOfMaps, mTermMaps.end() )
        {}
            
        /**
         * Set the number of equal-term-detection-maps this manager manages to the given value.
         * @param _newSize The value to set the number of equal-term-detection-maps this manager manages to.
         */
        void resize( size_t _newSize )
        {
            assert( !inUse() );
            mUsers( _newSize, nullptr );
            mTermMaps( _newSize, MapType() );
            mLeadingTerms( _newSize, mTermMaps.end() );
            mConstantTerms( _newSize, mTermMaps.end() );
        }
        
        /**
         * Reserve a free equal-term-detection-map.
         * @param _user The polynomial, which wants to reserve a equal-term-detection-map.
         * @param _expectedSize An upper bound of the expected number of terms, which are going to be added.
         * @return The id of the reserved equal-term-detection-map.
         */
        size_t getTermMapId( const Polynomial& _user, size_t _expectedSize )
        {
            std::lock_guard<std::mutex> lock(this->freshVarMutex);
            assert( mUsers.at( mNextMapId ) == nullptr );
            assert( mTermMaps.at( mNextMapId ).empty());
            assert( mLeadingTerms.at( mNextMapId ) == mTermMaps.end() );
            assert( mConstantTerms.at( mNextMapId ) == mTermMaps.end() );
            mUsers[mNextMapId] = &_user;
            mTermMaps[mNextMapId].reserve( _expectedSize );
            size_t result = mNextMapId;
            mNextMapId = (mNextMapId + 1) % mTermMaps.size();
            return result;
        }
        
        /**
         * Add the given term to the equal-term-detection-map with the given id.
         * @param _user The polynomial which uses this map.
         * @param _id The id of the map to add the term to.
         * @param _term The term to add.
         */
        void addTerm( const Polynomial& _user, size_t _id, std::shared_ptr<const Term<typename Polynomial::CoeffType>> _term )
        {
            assert( mUsers.at( _id ) == &_user );
            assert( mTermMaps[_id].size() < mTermMaps[_id].capacity() );
            auto res = mTermMaps[_id].emplace( _term->monomial(), _term );
            if( res.second )
            {
                if( _term->monomial() == nullptr )
                {
                    assert( mConstantTerms.at( _id ) == mTermMaps.end() );
                    mConstantTerms[_id] = res.first;
                }
                else if( mLeadingTerms[_id] == mTermMaps.end() || *(mLeadingTerms[_id]->first) < *_term->monomial() )
                {
                    mLeadingTerms[_id] = res.first;
                }
            }
            else
            {
                typename Polynomial::CoeffType tmp = res.first->second->coeff() + _term->coeff();
                if( tmp == 0 )
                {
                    mTermMaps[_id].erase( res.first );
                }
                else
                {
                    res.first->second = std::make_shared<const Term<typename Polynomial::CoeffType>>( tmp, _term->monomial() );
                }
            }
        }
        
        /**
         * Copy the terms in the equal-term-detection-map to the given term vector, while putting the leading term at the end of this vector
         * and the constant part at the beginning.
         * @sideeffect The given vector to store the terms into will be cleared an resized before adding terms to it.
         * @sideeffect The given id is released and the corresponding equal-term-detection-map is cleared.
         * @param _user The polynomial which uses this map.
         * @param _id The id of the map to add the term to.
         * @param _terms The vector of terms containing all terms of equal-term-detection-map with the given id.
         */
        void readTerms( const Polynomial& _user, size_t _id, std::vector<std::shared_ptr<const Term<typename Polynomial::CoeffType>>>& _terms )
        {
            assert( mUsers.at( _id ) == &_user );
            MapType& mcMap = mTermMaps[_id];
            _terms.clear();
            _terms.reserve( mcMap.size() );
            if( mLeadingTerms[_id] != mTermMaps.end() )
            {
                _terms.push_back( mLeadingTerms[_id]->second );
                mTermMaps.erase( mLeadingTerms[_id] );
                mLeadingTerms[_id] = mTermMaps.end();
            }
            std::shared_ptr<const Term<typename Polynomial::CoeffType>> constantPart = nullptr;
            if( mConstantTerms[_id] != mTermMaps.end() )
            {
                constantPart = mConstantTerms[_id]->second;
                mTermMaps.erase( mConstantTerms[_id] );
                mConstantTerms[_id] = mTermMaps.end();
            }
            for( auto iter = mTermMaps.begin(); iter != mTermMaps.end(); ++iter )
            {
                _terms.push_back( iter->second );
            }
            mTermMaps.clear();
            if( constantPart != nullptr )
            {
                _terms.push_back( constantPart );
            }
            mUsers[_id] = nullptr;
        }
        
     private:
         
        /**
         * @return true, if any equal-term-detection-map is in use.
         */
        bool inUse() const
        {
            for( int i = 0; i < mTermMaps.size(); ++i )
            {
                if( mUsers.at( i ) != nullptr ) return true;
                if( !mTermMaps.at( i ).empty() ) return true;
                if( mLeadingTerms.at( i ) != mTermMaps.end() ) return true;
                if( mConstantTerms.at( i ) != mTermMaps.end() ) return true;
            }
            return false;
        }
};

}

