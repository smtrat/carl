/**
 * Constraint.h
 * @author Florian Corzilius<corzilius@cs.rwth-aachen.de>
 * @author Sebastian Junges
 * @author Ulrich Loup
 * @since 2010-04-26
 * @version 2014-10-30
 */

#pragma once

//#define NDEBUG
#define CONSTRAINT_WITH_FACTORIZATION

#include "../config.h"
#include "../core/Relation.h"
#include "../core/Variables.h"
#include "../core/VariableInformation.h"
#include "../core/VariablesInformation.h"
#include "../core/polynomialfunctions/Complexity.h"
#include "../core/polynomialfunctions/Definiteness.h"
#include "../interval/Interval.h"
#include "../interval/IntervalEvaluation.h"
#include "../util/Common.h"
#include "config.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include <sstream>


namespace carl
{
    // Forward definition.
    template<typename Pol>
    class Constraint;
    
    template<typename Poly>
    using Constraints = std::set<Constraint<Poly>, carl::less<Constraint<Poly>, false>>;
    
    // Forward declaration
    template<typename Pol>
    class ConstraintPool;
    
    template<typename Pol>
    using VarInfo = VariableInformation<true, Pol>;
    
    template<typename Pol>
    using VarInfoMap = std::map<Variable, VarInfo<Pol>>;
    
    template<typename Pol, EnableIf<needs_cache<Pol>> = dummy>
    Pol makePolynomial( typename Pol::PolyType&& _poly );

    template<typename Pol, EnableIf<needs_cache<Pol>> = dummy>
    Pol makePolynomial( carl::Variable::Arg _var );
    
    template<typename Pol, EnableIf<needs_cache<Pol>> = dummy>
    Pol makePolynomial( const typename Pol::PolyType& _poly )
    {
        return makePolynomial<Pol>(typename Pol::PolyType(_poly));
    }

    template<typename Pol, DisableIf<needs_cache<Pol>> = dummy>
    Pol makePolynomial( carl::Variable::Arg _var )
    {
        return Pol( _var );
    }

    template<typename Pol, DisableIf<needs_cache<Pol>> = dummy>
    Pol makePolynomial( const typename Pol::PolyType& _poly )
    {
        return _poly;
    }
        
    /**
     * Represent a polynomial (in)equality against zero.
     */
    template<typename Pol>
    class ConstraintContent : public boost::intrusive::unordered_set_base_hook<>
    {
            friend class Constraint<Pol>;
            friend class ConstraintPool<Pol>;
        
        private:
            /// A unique id.
            /* const */ std::size_t mID;
            /// The relation symbol comparing the polynomial considered by this constraint to zero.
            const Relation mRelation;
            /// The polynomial which is compared by this constraint to zero.
            const Pol mLhs;
			/// The hash value.
            const carlVariables mVariables;
            /// Definiteness of the polynomial in this constraint.
            const Definiteness mLhsDefiniteness;
            /// Consistency
            const unsigned mConsistency;
            /// The factorization of the polynomial considered by this constraint.
            const std::size_t mHash;
            /// Weak pointer to itself
            std::weak_ptr<ConstraintContent<Pol>> mWeakPtr;
            /// A container which includes all variables occurring in the polynomial considered by this constraint.
            Factors<Pol> mFactorization;
            /// A map which stores information about properties of the variables in this constraint.
            VarInfoMap<Pol> mVarInfoMap;
            /// Mutex for access to variable information map.
            std::mutex mVarInfoMapMutex;
            /// Mutex for access to the factorization.
            std::mutex mFactorizationMutex;

            ConstraintContent() = delete;
            ConstraintContent(const ConstraintContent<Pol>&) = delete;
            
            template<typename P = Pol, EnableIf<needs_cache<P>> = dummy>
            ConstraintContent( std::size_t _id, typename Pol::PolyType&& _lhs, Relation _rel,  carl::carlVariables&& _vars, carl::Definiteness _definiteness, unsigned _consistent):
                ConstraintContent<Pol>::ConstraintContent( _id, std::move( makePolynomial<Pol>( std::move( _lhs ) ) ), _rel, std::move(_vars), _definiteness, _consistent )
            {}
            
            ConstraintContent( std::size_t _id, Pol&& _lhs, Relation _rel, carl::carlVariables&& _vars, carl::Definiteness _definiteness, unsigned _consistent );
                        
            /**
             * Initializes the stored factorization.
             */
            void initFactorization();
            
            void initVariableInformations() {
                VariablesInformation<false,Pol> varinfos = mLhs.template getVarInfo<false>();
				for (const auto& vi: varinfos) {
					mVarInfoMap.emplace_hint(mVarInfoMap.end(), vi.first, vi.second);
				}
            }
                           
        public:

            /**
             * Destructor.
             */
            ~ConstraintContent() noexcept {
                ConstraintPool<Pol>::getInstance().free(this);
            }

            /**
             * @return A hash value for this constraint.
             */
            std::size_t hash() const {
                return mHash;
            }
			
			std::size_t id() const {
				return mID;
			}
			Relation relation() const {
				return mRelation;
			}
			const auto& lhs() const {
				return mLhs;
			}

            unsigned isConsistent() const {
                return mConsistency;
            }
            
            /**
             * @param _variable The variable for which to determine the maximal degree.
             * @return The maximal degree of the given variable in this constraint content. (Monomial-wise)
             */
            uint maxDegree( const Variable& _variable ) const {
                auto varInfo = mVarInfoMap.find(_variable);
                if (varInfo == mVarInfoMap.end()) return 0;
                return varInfo->second.maxDegree();
            }
            
            /**
             * @return The maximal degree of all variables in this constraint. (Monomial-wise)
             */
            uint maxDegree() const {
                uint result = 0;
                for (const auto& var: mVariables) {
                    uint deg = maxDegree(var);
                    if (deg > result) result = deg;
                }
                return result;
            }
                  
    };
	
	/**
	 * @param lhs Left ConstraintContent
	 * @param rhs Right ConstraintContent
	 * @return `lhs == rhs`
	 */
	template<typename Pol>
	bool operator==(const ConstraintContent<Pol>& lhs, const ConstraintContent<Pol>& rhs) {
        return lhs.id() == rhs.id();
	}

    /**
     * Prints the representation of the given constraints on the given stream.
     * @param os The stream to print on.
     * @param cc The constraint to print.
     * @return The given stream after printing.
     */
    template<typename P>
    std::ostream& operator<<(std::ostream& os, const ConstraintContent<P>& cc) {
		return os << cc.lhs() << " " << cc.relation() << " 0";
    }

    /**
     * Represent a polynomial (in)equality against zero. Such an (in)equality
     * can be seen as an atomic formula/atom for the theory of real arithmetic.
     * Actually, this is just a (possibly) thread-safe wrapper with convenient
     * functions around the "ConstraintContent" class.
     */
    template<typename Pol>
    class Constraint
    {
            friend class ConstraintPool<Pol>;
            
        private:
            
            // Members.

            /// The real (in)equality is stored here.
            std::shared_ptr<ConstraintContent<Pol>> mpContent;
            
            explicit Constraint( std::shared_ptr<ConstraintContent<Pol>> _content );
            
            #ifdef THREAD_SAFE
            #define VARINFOMAP_LOCK_GUARD std::lock_guard<std::mutex> lock1( mpContent->mVarInfoMapMutex );
            #define FACTORIZATION_LOCK_GUARD std::lock_guard<std::mutex> lock1( mpContent->mFactorizationMutex );
            #define FACTORIZATION_LOCK mpContent->mFactorizationMutex.lock();
            #define FACTORIZATION_UNLOCK mpContent->mFactorizationMutex.unlock();
            #else
            #define VARINFOMAP_LOCK_GUARD
            #define FACTORIZATION_LOCK_GUARD
            #define FACTORIZATION_LOCK
            #define FACTORIZATION_UNLOCK
            #endif
            
        public:
            
            explicit Constraint( bool _valid = true );
            
            explicit Constraint( carl::Variable::Arg _var, Relation _rel, const typename Pol::NumberType& _bound = constant_zero<typename Pol::NumberType>::get() );
            
            explicit Constraint( const Pol& _lhs, Relation _rel );
            
            template<typename P = Pol, EnableIf<needs_cache<P>> = dummy>
            explicit Constraint( const typename P::PolyType& _lhs, Relation _rel );
            
            Constraint( const Constraint& _constraint );
            
            Constraint( Constraint&& _constraint ) noexcept;
                        
            Constraint& operator=( const Constraint& _constraint );
            
            Constraint& operator=( Constraint&& _constraint ) noexcept;
            
            /**
             * @return The considered polynomial being the left-hand side of this constraint.
             *          Hence, the right-hand side of any constraint is always 0.
             */
            const Pol& lhs() const
            {
                return mpContent->lhs();
            }

            /**
             * @return A container containing all variables occurring in the polynomial of this constraint.
             */
            const auto& variables() const
            {
                return mpContent->mVariables;
            }

			void gatherVariables(carlVariables& vars) const {
				vars.add(mpContent->mVariables.begin(), mpContent->mVariables.end());
			}

            /**
             * @return The relation symbol of this constraint.
             */
            Relation relation() const
            {
                return mpContent->relation();
            }

            /**
             * @return The unique id of this constraint.
             */
            size_t id() const
            {
                return mpContent->id();
            }

            /**
             * @return A hash value for this constraint.
             */
            size_t getHash() const
            {
                return mpContent->hash();
            }

            /**
             * @return true, if the polynomial p compared by this constraint has a proper factorization (!=p);
             *          false, otherwise.
             */
            bool hasFactorization() const
            {
                FACTORIZATION_LOCK
                if( mpContent->mFactorization.empty() )
                    mpContent->initFactorization();
                FACTORIZATION_UNLOCK
                return (mpContent->mFactorization.size() > 1);
            }

            /**
             * @return The factorization of the polynomial compared by this constraint.
             */
            const Factors<Pol>& factorization() const
            {
                FACTORIZATION_LOCK
                if( mpContent->mFactorization.empty() )
                    mpContent->initFactorization();
                FACTORIZATION_UNLOCK
                return mpContent->mFactorization;
            }

            /**
             * @return The constant part of the polynomial compared by this constraint.
             */
            typename Pol::NumberType constantPart() const
            {
                return mpContent->mLhs.constantPart();
            }
            
            /**
             * @param _variable The variable for which to determine the maximal degree.
             * @return The maximal degree of the given variable in this constraint. (Monomial-wise)
             */
            uint maxDegree( const Variable& _variable ) const
            {   
                VARINFOMAP_LOCK_GUARD
                return mpContent->maxDegree( _variable );
            }

            /**
             * @return The maximal degree of all variables in this constraint. (Monomial-wise)
             */
            uint maxDegree() const
            {
                VARINFOMAP_LOCK_GUARD
                return mpContent->maxDegree();
            }

            /**
             * @param _variable The variable for which to determine the minimal degree.
             * @return The minimal degree of the given variable in this constraint. (Monomial-wise)
             */
            uint minDegree( const Variable& _variable ) const
            {
                VARINFOMAP_LOCK_GUARD
                auto varInfo = mpContent->mVarInfoMap.find( _variable );
                if( varInfo == mpContent->mVarInfoMap.end() ) return 0;
                return varInfo->second.minDegree();
            }
            
            /**
             * @param _variable The variable for which to determine the number of occurrences.
             * @return The number of occurrences of the given variable in this constraint. (In 
             *          how many monomials of the left-hand side does the given variable occur?)
             */
            uint occurences( const Variable& _variable ) const
            {
                VARINFOMAP_LOCK_GUARD
                auto varInfo = mpContent->mVarInfoMap.find( _variable );
                if( varInfo == mpContent->mVarInfoMap.end() ) return 0;
                return varInfo->second.occurence();
            }
            
            /**
             * @param _variable The variable to find variable information for.
			 * @param _withCoefficients
             * @return The whole variable information object.
             * Note, that if the given variable is not in this constraints, this method fails.
             * Furthermore, the variable information returned do provide coefficients only, if
             * the given flag _withCoefficients is set to true.
             */
            const VarInfo<Pol>& varInfo( const Variable& _variable, bool _withCoefficients = false ) const
            {
                VARINFOMAP_LOCK_GUARD
                auto varInfo = mpContent->mVarInfoMap.find( _variable );
                assert( varInfo != mpContent->mVarInfoMap.end() );
                if( _withCoefficients && !varInfo->second.hasCoeff() )
                {
                    varInfo->second = mpContent->mLhs.template getVarInfo<true>( _variable );
                }
                return varInfo->second;
            }
			
			bool relationIsStrict() const {
                return isStrict(mpContent->mRelation);
			}
			bool relationIsWeak() const {
				return isWeak(mpContent->mRelation);
			}
            
            /**
             * Checks if the given variable occurs in the constraint.
             * @param _var  The variable to check for.
             * @return true, if the given variable occurs in the constraint;
             *          false, otherwise.
             */
            bool hasVariable( const Variable& _var ) const
            {
				return mpContent->mVariables.has(_var);
            }
            
            /**
             * @return true, if it contains only integer valued variables.
             */
            bool integerValued() const
            {
				return mpContent->mVariables.size() == mpContent->mVariables.integer().size();
            }
            
            /**
             * @return true, if it contains only real valued variables.
             */
            bool realValued() const
            {
				return mpContent->mVariables.size() == mpContent->mVariables.real().size();
            }
            
            /**
             * Checks if this constraints contains an integer valued variable.
             * @return true, if it does;
             *          false, otherwise.
             */
            bool hasIntegerValuedVariable() const
            {
                return !mpContent->mVariables.integer().empty();
            }
            
            /**
             * Checks if this constraints contains an real valued variable.
             * @return true, if it does;
             *          false, otherwise.
             */
            bool hasRealValuedVariable() const
            {
                return !mpContent->mVariables.real().empty();
            }
            
            /**
             * @return true, if this constraint is a bound.
             */
            bool isBound(bool negated = false) const
            {
				if (mpContent->mVariables.size() != 1 || maxDegree(mpContent->mVariables.as_vector()[0]) != 1) return false;
				if (negated) {
					return mpContent->mRelation != Relation::EQ;
				} else {
					return mpContent->mRelation != Relation::NEQ;
				}
            }
            
            /**
             * @return true, if this constraint is a lower bound.
             */
            bool isLowerBound() const
            {
                if( isBound() )
                {
                    if( mpContent->mRelation == Relation::EQ ) return true;
                    const typename Pol::NumberType& coeff = mpContent->mLhs.lterm().coeff();
                    if( coeff < 0 )
                        return (mpContent->mRelation == Relation::LEQ || mpContent->mRelation == Relation::LESS );
                    else
                    {
                        assert( coeff > 0 );
                        return (mpContent->mRelation == Relation::GEQ || mpContent->mRelation == Relation::GREATER );
                    }
                }
                return false;
            }
            
            /**
             * @return true, if this constraint is an upper bound.
             */
            bool isUpperBound() const
            {
                if( isBound() )
                {
                    if( mpContent->mRelation == Relation::EQ ) return true;
                    const typename Pol::NumberType& coeff = mpContent->mLhs.lterm().coeff();
                    if( coeff > 0 )
                        return (mpContent->mRelation == Relation::LEQ || mpContent->mRelation == Relation::LESS );
                    else
                    {
                        assert( coeff < 0 );
                        return (mpContent->mRelation == Relation::GEQ || mpContent->mRelation == Relation::GREATER );
                    }
                }
                return false;
            }
            
            /**
             * @return An approximation of the complexity of this constraint.
             */
            size_t complexity() const
            {
				return 1 + carl::complexity(mpContent->mLhs);
            }
            
            /**
             * Checks whether the given assignment satisfies this constraint.
             * @param _assignment The assignment.
             * @return 1, if the given assignment satisfies this constraint.
             *          0, if the given assignment contradicts this constraint.
             *          2, otherwise (possibly not defined for all variables in the constraint,
             *                       even then it could be possible to obtain the first two results.)
             */
            unsigned satisfiedBy( const EvaluationMap<typename Pol::NumberType>& _assignment ) const;
            
            /**
             * Checks, whether the constraint is consistent.
             * It differs between, containing variables, consistent, and inconsistent.
             * @return 0, if the constraint is not consistent.
             *          1, if the constraint is consistent.
             *          2, if the constraint still contains variables.
             */
            unsigned isConsistent() const
            {
                return mpContent->isConsistent();
            }
            
            /**
             * Checks whether this constraint is consistent with the given assignment from 
             * the its variables to interval domains.
             * @param _solutionInterval The interval domains of the variables.
             * @return 1, if this constraint is consistent with the given intervals;
             *          0, if this constraint is not consistent with the given intervals;
             *          2, if it cannot be decided whether this constraint is consistent with the given intervals.
             */
            unsigned consistentWith( const EvaluationMap<Interval<double>>& _solutionInterval ) const;
            
            /**
             * Checks whether this constraint is consistent with the given assignment from 
             * the its variables to interval domains.
             * @param _solutionInterval The interval domains of the variables.
             * @param _stricterRelation This relation is set to a relation R such that this constraint and the given variable bounds
             *                           imply the constraint formed by R, comparing this constraint's left-hand side to zero.
             * @return 1, if this constraint is consistent with the given intervals;
             *          0, if this constraint is not consistent with the given intervals;
             *          2, if it cannot be decided whether this constraint is consistent with the given intervals.
             */
            unsigned consistentWith( const EvaluationMap<Interval<double>>& _solutionInterval, Relation& _stricterRelation ) const;

			/**
			 * Checks whether the given interval assignment may fulfill the constraint.
			 * Note that the assignment must be complete.
			 * There are three possible outcomes:
			 * - True (4), i.e. all actual assignments satisfy the constraint: $p ~_\alpha 0 \Leftrightarrow True$
			 * - Maybe (3), i.e. some actual assignments satisfy the constraint: $p ~_\alpha 0 \Leftrightarrow ?$
			 * - Not null (2), i.e. all assignments that make the constraint evaluate not to zero satisfy the constraint: $p ~_\alpha 0 \Leftrightarrow p \neq 0$
			 * - Null (1), i.e. only assignments that make the constraint evaluate to zero satisfy the constraint: $p ~_\alpha 0 \Leftrightarrow p_\alpha = 0$
			 * - False (0), i.e. no actual assignment satisfies the constraint: $p ~_\alpha 0 \Leftrightarrow False$
			 * @param _assignment Variable assignment.
			 * @return 0, 1 or 2.
			 */
			unsigned evaluate(const EvaluationMap<Interval<typename carl::UnderlyingNumberType<Pol>::type>>& _assignment) const;

            /**
             * @param _var The variable to check the size of its solution set for.
             * @return true, if it is easy to decide whether this constraint has a finite solution set
             *                in the given variable;
             *          false, otherwise.
             */
            bool hasFinitelyManySolutionsIn( const Variable& _var ) const;
                        
            /**
             * Calculates the coefficient of the given variable with the given degree. Note, that it only
             * computes the coefficient once and stores the result.
             * @param _var The variable for which to calculate the coefficient.
             * @param _degree The according degree of the variable for which to calculate the coefficient.
             * @return The ith coefficient of the given variable, where i is the given degree.
             */
            Pol coefficient( const Variable& _var, uint _degree ) const;
            
            Constraint negation() const {
                CARL_LOG_DEBUG("carl.formula.constraint", "negation of " << *this << " is " << Constraint(lhs(), carl::inverse(relation())));
                return Constraint(lhs(), carl::inverse(relation()));
            }

            /**
             * If this constraint represents a substitution (equation, where at least one variable occurs only linearly),
             * this method detects a (there could be various possibilities) corresponding substitution variable and term.
             * @param _substitutionVariable Is set to the substitution variable, if this constraint represents a substitution.
             * @param _substitutionTerm Is set to the substitution term, if this constraint represents a substitution.
             * @return true, if this constraints represents a substitution;
             *         false, otherwise.
             */
            bool getSubstitution( Variable& _substitutionVariable, Pol& _substitutionTerm, bool _negated = false, const Variable& _exclude = carl::Variable::NO_VARIABLE ) const;
			
			bool getAssignment(Variable& _substitutionVariable, typename Pol::NumberType& _substitutionValue) const;

            /**
             * Determines whether the constraint is pseudo-boolean.
             *
             * @return True if this constraint is pseudo-boolean. False otherwise.
             */
            bool isPseudoBoolean() const;

            /**
             * Prints the properties of this constraints on the given stream.
             * @param _out The stream to print on.
             */
            void printProperties( std::ostream& _out = std::cout ) const;
		
		template<typename P>
		friend bool operator==(const Constraint<P>& lhs, const Constraint<P>& rhs);
		template<typename P>
		friend bool operator!=(const Constraint<P>& lhs, const Constraint<P>& rhs);
    };
	
	/**
	 * @param lhs Left constraint
	 * @param rhs Right constraint
	 * @return `lhs == rhs`
	 */
	template<typename P>
	bool operator==(const Constraint<P>& lhs, const Constraint<P>& rhs) {
		return lhs.mpContent == rhs.mpContent;
	}
	
	/**
	 * @param lhs Left constraint
	 * @param rhs Right constraint
	 * @return `lhs != rhs`
	 */
	template<typename P>
	bool operator!=(const Constraint<P>& lhs, const Constraint<P>& rhs) {
		return lhs.mpContent != rhs.mpContent;
	}
	
	/**
	 * @param lhs Left constraint
	 * @param rhs Right constraint
	 * @return `lhs < rhs`
	 */
	template<typename P>
	bool operator<(const Constraint<P>& lhs, const Constraint<P>& rhs) {
		assert(lhs.id() != 0 && rhs.id() != 0);
		return lhs.id() < rhs.id();
	}
	
	/**
	 * @param lhs Left constraint
	 * @param rhs Right constraint
	 * @return `lhs > rhs`
	 */
	 template<typename P>
	 bool operator>(const Constraint<P>& lhs, const Constraint<P>& rhs) {
		assert(lhs.id() != 0 && rhs.id() != 0);
		return lhs.id() > rhs.id();
	 }
	
	 /**
 	 * @param lhs Left constraint
 	 * @param rhs Right constraint
 	 * @return `lhs <= rhs`
 	 */
	 template<typename P>
	 bool operator<=(const Constraint<P>& lhs, const Constraint<P>& rhs) {
		assert(lhs.id() != 0 && rhs.id() != 0);
		return lhs.id() <= rhs.id();
	 }
	
	 /**
 	 * @param lhs Left constraint
 	 * @param rhs Right constraint
 	 * @return `lhs >= rhs`
 	 */
	 template<typename P>
	 bool operator>=(const Constraint<P>& lhs, const Constraint<P>& rhs) {
		assert(lhs.id() != 0 && rhs.id() != 0);
		return lhs.id() >= rhs.id();
	 }

    
    const signed A_IFF_B = 2;
    const signed A_IMPLIES_B = 1;
    const signed B_IMPLIES_A = -1;
    const signed NOT__A_AND_B = -2;
    const signed A_AND_B__IFF_C = -3;
    const signed A_XOR_B = -4;
            
    /**
     * Compares _constraintA with _constraintB.
     * @return  2, if it is easy to decide that _constraintA and _constraintB have the same solutions. _constraintA = _constraintB
     *           1, if it is easy to decide that _constraintB includes all solutions of _constraintA;   _constraintA -> _constraintB
     *          -1, if it is easy to decide that _constraintA includes all solutions of _constraintB;   _constraintB -> _constraintA
     *          -2, if it is easy to decide that _constraintA has no solution common with _constraintB; not(_constraintA and _constraintB)
     *          -3, if it is easy to decide that _constraintA and _constraintB can be intersected;      _constraintA and _constraintB = _constraintC
     *          -4, if it is easy to decide that _constraintA is the inverse of _constraintB;           _constraintA xor _constraintB
     *           0, otherwise.
     */
    template<typename Pol>
    signed compare( const Constraint<Pol>& _constraintA, const Constraint<Pol>& _constraintB )
    {
        /*
         * Check whether it holds that 
         * 
         *                      _constraintA  =  a_1*m_1+...+a_k*m_k + c ~ 0
         * and 
         *                      _constraintB  =  b_1*m_1+...+b_k*m_k + d ~ 0, 
         * 
         * where a_1,..., a_k, b_1,..., b_k, c, d are rational coefficients, 
         *       m_1,..., m_k are non-constant monomials and 
         *       exists a rational g such that 
         * 
         *                   a_i = g * b_i for all 1<=i<=k 
         *              or   b_i = g * a_i for all 1<=i<=k 
         */
        typename Pol::NumberType one_divided_by_a = _constraintA.lhs().coprimeFactorWithoutConstant();
        typename Pol::NumberType one_divided_by_b = _constraintB.lhs().coprimeFactorWithoutConstant();
        typename Pol::NumberType c = _constraintA.lhs().constantPart();
        typename Pol::NumberType d = _constraintB.lhs().constantPart();
        assert( carl::isOne(carl::getNum(carl::abs(one_divided_by_b))) );
        Pol tmpA = (_constraintA.lhs() - c) * one_divided_by_a;
        Pol tmpB = (_constraintB.lhs() - d) * one_divided_by_b;
//        std::cout << "tmpA = " << tmpA << std::endl;
//        std::cout << "tmpB = " << tmpB << std::endl;
        if( tmpA != tmpB ) return 0;
        bool termACoeffGreater = false;
        bool signsDiffer = (one_divided_by_a < carl::constant_zero<typename Pol::NumberType>::get()) != (one_divided_by_b < carl::constant_zero<typename Pol::NumberType>::get());
        typename Pol::NumberType g;
        if( carl::getDenom(one_divided_by_a) > carl::getDenom(one_divided_by_b) )
        {
            g = typename Pol::NumberType(carl::getDenom(one_divided_by_a))/carl::getDenom(one_divided_by_b);
            if( signsDiffer )
                g = -g;
            termACoeffGreater = true;
            d *= g;
        }
        else
        {
            g = typename Pol::NumberType(carl::getDenom(one_divided_by_b))/carl::getDenom(one_divided_by_a);
            if( signsDiffer )
                g = -g;
            c *= g;
        }
        // Apply the multiplication by a negative g to the according relation symbol, which
        // has to be turned around then.
        Relation relA = _constraintA.relation();
        Relation relB = _constraintB.relation();
        if( g < 0 )
        {
            if( termACoeffGreater )
            {
                switch( relB )
                {
                    case Relation::LEQ:
                        relB = Relation::GEQ;
                        break;
                    case Relation::GEQ:
                        relB = Relation::LEQ;
                        break;
                    case Relation::LESS:
                        relB = Relation::GREATER;
                        break;
                    case Relation::GREATER:
                        relB = Relation::LESS;
                        break;
                    default:
                        break;
                }
            }
            else
            {
                switch( relA )
                {
                    case Relation::LEQ:
                        relA = Relation::GEQ;
                        break;
                    case Relation::GEQ:
                        relA = Relation::LEQ;
                        break;
                    case Relation::LESS:
                        relA = Relation::GREATER;
                        break;
                    case Relation::GREATER:
                        relA = Relation::LESS;
                        break;
                    default:
                        break;
                }   
            }
        }
//        std::cout << "c = " << c << std::endl;
//        std::cout << "d = " << d << std::endl;
//        std::cout << "g = " << g << std::endl;
//        std::cout << "relA = " << relA << std::endl;
//        std::cout << "relB = " << relB << std::endl;
        // Compare the adapted constant parts.
        switch( relB )
        {
            case Relation::EQ:
                switch( relA )
                {
                    case Relation::EQ: // p+c=0  and  p+d=0
                        if( c == d ) return A_IFF_B; 
                        return NOT__A_AND_B;
                    case Relation::NEQ: // p+c!=0  and  p+d=0
                        if( c == d ) return A_XOR_B;
                        return B_IMPLIES_A;
                    case Relation::LESS: // p+c<0  and  p+d=0
                        if( c < d ) return B_IMPLIES_A;
                        return NOT__A_AND_B;
                    case Relation::GREATER: // p+c>0  and  p+d=0
                        if( c > d ) return B_IMPLIES_A;
                        return NOT__A_AND_B;
                    case Relation::LEQ: // p+c<=0  and  p+d=0
                        if( c <= d ) return B_IMPLIES_A;
                        return NOT__A_AND_B;
                    case Relation::GEQ: // p+c>=0  and  p+d=0
                        if( c >= d ) return B_IMPLIES_A;
                        return NOT__A_AND_B;
                    default:
                        return false;
                }
            case Relation::NEQ:
                switch( relA )
                {
                    case Relation::EQ: // p+c=0  and  p+d!=0
                        if( c == d ) return A_XOR_B;
                        return A_IMPLIES_B;
                    case Relation::NEQ: // p+c!=0  and  p+d!=0
                        if( c == d ) return A_IFF_B;
                        return 0;
                    case Relation::LESS: // p+c<0  and  p+d!=0
                        if( c >= d ) return A_IMPLIES_B;
                        return 0;
                    case Relation::GREATER: // p+c>0  and  p+d!=0
                        if( c <= d ) return A_IMPLIES_B;
                        return 0;
                    case Relation::LEQ: // p+c<=0  and  p+d!=0
                        if( c > d ) return A_IMPLIES_B;
                        if( c == d ) return A_AND_B__IFF_C;
                        return 0;
                    case Relation::GEQ: // p+c>=0  and  p+d!=0
                        if( c < d ) return A_IMPLIES_B;
                        if( c == d ) return A_AND_B__IFF_C;
                        return 0;
                    default:
                        return 0;
                }
            case Relation::LESS:
                switch( relA )
                {
                    case Relation::EQ: // p+c=0  and  p+d<0
                        if( c > d ) return A_IMPLIES_B;
                        return NOT__A_AND_B;
                    case Relation::NEQ: // p+c!=0  and  p+d<0
                        if( c <= d ) return B_IMPLIES_A;
                        return 0;
                    case Relation::LESS: // p+c<0  and  p+d<0
                        if( c == d ) return A_IFF_B;
                        if( c < d ) return B_IMPLIES_A;
                        return A_IMPLIES_B;
                    case Relation::GREATER: // p+c>0  and  p+d<0
                        if( c <= d ) return NOT__A_AND_B;
                        return 0;
                    case Relation::LEQ: // p+c<=0  and  p+d<0
                        if( c > d ) return A_IMPLIES_B;
                        return B_IMPLIES_A;
                    case Relation::GEQ: // p+c>=0  and  p+d<0
                        if( c < d ) return NOT__A_AND_B;
                        if( c == d ) return A_XOR_B;
                        return 0;
                    default:
                        return 0;
                }
            case Relation::GREATER:
            {
                switch( relA )
                {
                    case Relation::EQ: // p+c=0  and  p+d>0
                        if( c < d ) return A_IMPLIES_B;
                        return NOT__A_AND_B;
                    case Relation::NEQ: // p+c!=0  and  p+d>0
                        if( c >= d ) return B_IMPLIES_A;
                        return 0;
                    case Relation::LESS: // p+c<0  and  p+d>0
                        if( c >= d ) return NOT__A_AND_B;
                        return 0;
                    case Relation::GREATER: // p+c>0  and  p+d>0
                        if( c == d ) return A_IFF_B;
                        if( c > d ) return B_IMPLIES_A;
                        return A_IMPLIES_B;
                    case Relation::LEQ: // p+c<=0  and  p+d>0
                        if( c > d ) return NOT__A_AND_B;
                        if( c == d ) return A_XOR_B;
                        return 0;
                    case Relation::GEQ: // p+c>=0  and  p+d>0
                        if( c > d ) return B_IMPLIES_A;
                        return A_IMPLIES_B;
                    default:
                        return 0;
                }
            }
            case Relation::LEQ:
            {
                switch( relA )
                {
                    case Relation::EQ: // p+c=0  and  p+d<=0
                        if( c >= d ) return A_IMPLIES_B;
                        return NOT__A_AND_B;
                    case Relation::NEQ: // p+c!=0  and  p+d<=0
                        if( c < d ) return B_IMPLIES_A;
                        if( c == d ) return A_AND_B__IFF_C;
                        return 0;
                    case Relation::LESS: // p+c<0  and  p+d<=0
                        if( c < d ) return B_IMPLIES_A;
                        return A_IMPLIES_B;
                    case Relation::GREATER: // p+c>0  and  p+d<=0
                        if( c < d ) return NOT__A_AND_B;
                        if( c == d ) return A_XOR_B;
                        return 0;
                    case Relation::LEQ: // p+c<=0  and  p+d<=0
                        if( c == d ) return A_IFF_B;
                        if( c < d ) return B_IMPLIES_A;
                        return A_IMPLIES_B;
                    case Relation::GEQ: // p+c>=0  and  p+d<=0
                        if( c < d ) return NOT__A_AND_B;
                        if( c == d ) return A_AND_B__IFF_C;
                        return 0;
                    default:
                        return 0;
                }
            }
            case Relation::GEQ:
            {
                switch( relA )
                {
                    case Relation::EQ: // p+c=0  and  p+d>=0
                        if( c <= d ) return A_IMPLIES_B;
                        return NOT__A_AND_B;
                    case Relation::NEQ: // p+c!=0  and  p+d>=0
                        if( c > d ) return B_IMPLIES_A;
                        if( c == d ) return A_AND_B__IFF_C;
                        return 0;
                    case Relation::LESS: // p+c<0  and  p+d>=0
                        if( c > d ) return NOT__A_AND_B;
                        if( c == d ) return A_XOR_B;
                        return 0;
                    case Relation::GREATER: // p+c>0  and  p+d>=0
                        if( c < d ) return B_IMPLIES_A;
                        return A_IMPLIES_B;
                    case Relation::LEQ: // p+c<=0  and  p+d>=0
                        if( c > d ) return NOT__A_AND_B;
                        if( c == d ) return A_AND_B__IFF_C;
                        return 0;
                    case Relation::GEQ: // p+c>=0  and  p+d>=0
                        if( c == d ) return A_IFF_B;
                        if( c < d ) return A_IMPLIES_B;
                        return B_IMPLIES_A;
                    default:
                        return 0;
                }
            }
            default:
                return 0;
        }
    }
    
    /**
     * Prints the given constraint on the given stream.
     * @param os The stream to print the given constraint on.
     * @param c The formula to print.
     * @return The stream after printing the given constraint on it.
     */
    template<typename Poly>
    std::ostream& operator<<(std::ostream& os, const Constraint<Poly>& c) {
		return os << c.lhs() << " " << c.relation() << " 0";
	}
    
}    // namespace carl

#define CONSTRAINT_HASH( _lhs, _rel, _type ) (size_t)((size_t)(std::hash<_type>()( _lhs ) << 3) ^ (size_t)_rel)

namespace std
{
    /**
     * Implements std::hash for constraint contents.
     */
    template<typename Pol>
    struct hash<carl::ConstraintContent<Pol>> {
        /**
         * @param _constraintContent The constraint content to get the hash for.
         * @return The hash of the given constraint content.
         */
        std::size_t operator()( const carl::ConstraintContent<Pol>& _constraintContent ) const 
        {
            return _constraintContent.hash();
        }
    };
    
    /**
     * Implements std::hash for constraints.
     */
    template<typename Pol>
    struct hash<carl::Constraint<Pol>> {
        /**
         * @param _constraint The constraint to get the hash for.
         * @return The hash of the given constraint.
         */
        std::size_t operator()( const carl::Constraint<Pol>& _constraint ) const 
        {
            return _constraint.getHash();
        }
    };
    
    /**
     * Implements std::hash for vectors of constraints.
     */
    template<typename Pol>
    struct hash<std::vector<carl::Constraint<Pol>>> {
        /**
         * @param _arg The vector of constraints to get the hash for.
         * @return The hash of the given vector of constraints.
         */
        std::size_t operator()( const std::vector<carl::Constraint<Pol>>& _arg ) const
        {
            std::size_t result = 0;
			for (const auto& c: _arg) {
                result <<= 5;
                result ^= c.id();
            }
            return result;
        }
    };
} // namespace std

#include "Constraint.tpp"
