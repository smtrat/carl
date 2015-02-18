/**
 * @file SortManager.h
 * @author Gereon Kremer <gereon.kremer@cs.rwth-aachen.de>
 * @author Florian Corzilius <corzilius@cs.rwth-aachen.de>
 * @since 2014-10-30
 * @version 2014-10-30
 */

#pragma once

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "../core/logging.h"
#include "../util/Singleton.h"
#include "../util/Common.h"
#include "../io/streamingOperators.h"
#include "Sort.h"

namespace carl
{

/**
 * Implements a manager for sorts, containing the actual contents of these sort and allocating their ids.
 */
class SortManager : public Singleton<SortManager>
{
	
		friend Singleton<SortManager>;

	public:
		
		/// The actual content of a sort.
		struct SortContent
		{
			/// The sort's name.
			std::string name;
			union
			{
				/// The sort's carl variable type, if this sort is interpreted.
				VariableType type;
				/// The sort's argument types, if this sort is not interpreted. It is nullptr, if the sort's arity is zero.
				std::vector<Sort>* parameters;
			};
			/// A flag specifying whether the sort is interpreted (e.g., the real or integer domain) or not.
			bool interpreted;
			/// The sort's indices. A sort can be indexed with the "_" operator. It is nullptr, if no indices are present.
			std::vector<std::size_t>* indices;
			
			SortContent() = delete; // The default constructor is disabled.
			
			/**
			 * Constructs a sort content.
			 * @param _name The name of the sort content to construct.
			 * @param _type The carl variable type of the sort content to construct.
			 */
			explicit SortContent( const std::string& _name, VariableType _type ):
				name( _name ),
				type( _type ),
				interpreted( true )
			{}
			
			/**
			 * Constructs a sort content.
			 * @param _name The name of the sort content to construct.
			 */
			explicit SortContent( const std::string& _name ):
				name( _name ),
				parameters( nullptr ),
				interpreted( false ),
				indices( nullptr )
			{}
			
			/**
			 * Constructs a sort content.
			 * @param _name The name  of the sort content to construct.
			 * @param _parameters The sorts of the arguments of the sort content to construct.
			 */
			explicit SortContent( const std::string& _name, std::vector<Sort>&& _parameters ):
				name( _name ),
				parameters( new std::vector<Sort>( std::move( _parameters ) ) ),
				interpreted( false ),
				indices( nullptr )
			{}
			
			SortContent(const SortContent& sc):
				name(sc.name),
				interpreted(sc.interpreted),
				indices(nullptr)
			{
				if (interpreted) type = sc.type;
				else if (sc.parameters != nullptr) parameters = new std::vector<Sort>(*sc.parameters);
				if (sc.indices != nullptr) indices = new std::vector<std::size_t>(*sc.indices);
			}
			
			/// Destructs a sort content.
			~SortContent() {
				if (!interpreted) delete parameters;
				delete indices;
			}
			
			/**
			 * @param _sc The sort content to compare with.
			 * @return true, if this sort content is less than the given one.
			 */
			bool operator<( const SortContent& _sc ) const
			{
				if (name != _sc.name) return name < _sc.name;
				if (interpreted != _sc.interpreted) return interpreted;
				if (interpreted) {
					if (type != _sc.type) return type < _sc.type;
				}
				if (!interpreted) {
					if (parameters == nullptr && _sc.parameters != nullptr) return true;
					if (parameters != nullptr && _sc.parameters == nullptr) return false;
					if (parameters != nullptr && _sc.parameters != nullptr) {
						if (*parameters != *_sc.parameters) return *parameters < *_sc.parameters;
					}
				}
				if (indices == nullptr && _sc.indices != nullptr) return true;
				if (indices != nullptr && _sc.indices == nullptr) return false;
				if (indices != nullptr && _sc.indices != nullptr) {
					if (*indices != *_sc.indices) return *indices < *_sc.indices;
				}
				return false;
			}
		};
		
		/// The type of a sort template, define by define-sort.
		typedef std::pair<std::vector<std::string>, Sort> SortTemplate;
		
	private:
		// Members.

		/// Stores all instantiated sorts and maps them to their unique id.
		PointerMap<SortContent, Sort::IDType> mSortcontentIdMap;
		/// Maps the unique ids to the sort content.
		std::vector<const SortContent*> mSorts;
		/// Stores all sort declarations invoked by a declare-sort.
		std::map<std::string, unsigned> mDeclarations;
		/// Stores all sort definitions invoked by a define-sort.
		std::map<std::string, SortTemplate> mDefinitions;
		/// Stores all interpreted sorts.
		std::map<VariableType, Sort> mInterpretedSorts;

		/**
		 * Constructs a sort manager.
		 */
		SortManager():
			mSortcontentIdMap(),
			mSorts(),
			mDeclarations(),
			mDefinitions()
		{
			mSorts.emplace_back( nullptr ); // default value
		}
		
		/**
		 * Tries to add the given sort content to the so far stored sort contents. If it has already been stored,
		 * if will be deleted and the id of the already existing sort content will be used to create the returned sort.
		 * @param _sc The sort content to store.
		 * @return The sort corresponding to the given sort content.
		 */
		Sort getSort( const SortContent* _sc );

	public:
		~SortManager() {
			for (auto& s: mSorts) delete s;
			mSorts.clear();
		}
		
		/**
		 * @param _sort A sort.
		 * @return The name if the given sort.
		 */
		const std::string& getName( const Sort& _sort ) const;
		
		const std::vector<std::size_t>* getIndices(const Sort& _sort) const {
			assert( _sort.id() != 0 );
			assert( _sort.id() < mSorts.size() );
			return mSorts[_sort.id()]->indices;
		}

		/**
		 * Prints the given sort on the given output stream.
		 * @param _out The output stream to print the given sort on.
		 * @param _sort The sort to print.
		 * @return The output stream after printing the given sort on it.
		 */
		std::ostream& print( std::ostream& _out, const Sort& _sort ) const;
		
		/**
		 * Instantiates the given sort according to the mapping of sort names to sorts as declared by the given map.
		 * @param _sort The sort to replace sorts by sorts in.
		 * @param _parameters The map of sort names to sorts.
		 * @return The resulting sort.
		 */
		Sort replace( const Sort& _sort, const std::map<std::string, Sort>& _parameters );

		/**
		 * Adds a sort declaration.
		 * @param _name The name of the declared sort.
		 * @param _arity The arity of the declared sort.
		 * @return true, if the given sort declaration has not been added before;
		 *		 false, otherwise.
		 */
		bool declare( const std::string& _name, unsigned _arity );

		/**
		 * Adds a sort template definitions.
		 * @param _name The name of the defined sort template.
		 * @param _params The template parameter of the defined sort.
		 * @param _sort The sort to instantiate into.
		 * @return true, if the given sort template definition has not been added before;
		 *		 false, otherwise.
		 */
		bool define( const std::string& _name, const std::vector<std::string>& _params, const Sort& _sort );

		/**
		 * @param _name The name of the sort declaration or sort template definition.
		 * @return The aritiy of the given sort declaration or sort template definition.
		 */
		std::size_t arity( const std::string& _name ) const;
		
		/**
		 * @param _sort The sort to get the arity for.
		 * @return The arity of the given sort.
		 */
		std::size_t getArity( const Sort& _sort ) const;

		/**
		 * Adds an interpreted sort.
		 * @param _name The name of the sort to add.
		 * @param type The carl variable type of the sort to add.
		 * @return The resulting sort.
		 */
		Sort interpretedSort( const std::string& _name, VariableType type );

		/**
		 * Retrieves an interpreted sort.
		 * Asserts, that the sort was already registered via interpretedSort(std::string, VariableType).
		 * @param vt VariableType.
		 * @return Sort corresponding to vt.
		 */
		const Sort& interpretedSort( VariableType vt ) const
		{
			auto it = mInterpretedSorts.find(vt);
			assert(it != mInterpretedSorts.end());
			return it->second;
		}

		/**
		 * @param _sort A sort.
		 * @return true, if the given sort is interpreted.
		 */
		bool isInterpreted( const Sort& _sort ) const
		{
			assert( _sort.id() != 0 );
			assert( _sort.id() < mSorts.size() );
			return mSorts[_sort.id()]->interpreted;
		}

		/**
		 * @param _sort A sort, which must be interpreted.
		 * @return The interpreted type/sort of the given sort.
		 */
		VariableType interpretedType( const Sort& _sort ) const
		{
			assert( _sort.id() != 0 );
			assert( _sort.id() < mSorts.size() );
			const SortContent& sc = *mSorts[_sort.id()];
			assert( sc.interpreted );
			return sc.type;
		}
		
		Sort instantiate(const std::string& name, const std::vector<Sort>& params);
		
		Sort index(const Sort& sort, const std::vector<std::size_t>& indices);

		/**
		 * Gets the sort with arity zero (thus it is maybe interpreted) corresponding the given name.
		 * @param _name The name of the sort to get.
		 * @return The resulting sort.
		 */
		Sort getSort(const std::string& _name);

		/**
		 * Gets the sort with arity greater than zero corresponding the given name and having the arguments
		 * of the given sorts.
		 * @param _name The name of the sort to get.
		 * @param _params The sort of the arguments of the sort to get.
		 * @return The resulting sort.
		 */
		Sort getSort(const std::string& _name, const std::vector<Sort>& _params);
		
		Sort getSort(const std::string& _name, const std::vector<std::size_t>& _indices);
		
		Sort getSort(const std::string& _name, const std::vector<std::size_t>& _indices, const std::vector<Sort>& _params);
};

/**
 * Gets the sort specified by the arguments.
 * Forwards to SortManager::newSort().
 */
template<typename... Args>
inline Sort getSort(Args&&... args) {
	return SortManager::getInstance().getSort(std::forward<Args>(args)...);
}

/*inline Sort getSort( const std::string& _name ) {
	return SortManager::getInstance().getSort(_name);
}
inline Sort getSort( const std::string& _name, const std::vector<Sort>& _params ) {
	return SortManager::getInstance().getSort(_name, _params);
}*/

}
