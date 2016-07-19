/* 
 * File:   ThomEncoding.h
 * Author: tobias
 *
 * Created on 24. April 2016, 20:04
 */

#pragma once


#include "../core/UnivariatePolynomial.h"
#include "../core/Sign.h"

#include "ThomUtil.h"
#include "TarskiQuery/TarskiQuery.h"


namespace carl {



/*
 * A class for both one- and multidimensional thom encodings.
 */
template<typename Coeff>
class ThomEncoding {
        
private:
                    
        /*
         * Stores a pointer to the polynomial
         */
        std::shared_ptr<MultivariatePolynomial<Coeff>> p; 
        
        /*
         * The list of sign conditions realized by the derivatives on the represented root
         * We only need to store the sign up to that of p^(deg(P) - 1) since p^(deg(P)) is constant
         */
        SignCondition signs;
        
        /*
         * the main variable on this level. The derivatives are partial w.r.t. this variable.
         */
        Variable mainVar;
        
        /*
         * The point of dimension d-1 represented by the encoding on the level below or nullptr if n.a.
         */
        std::shared_ptr<ThomEncoding<Coeff>> point;
        
public:
        
        //////////////////
        // Constructors //
        //////////////////
              
        /*
         * Constructs a univariate thom encoding with the given polynomial and the sign conditions
         */
        ThomEncoding(
                std::shared_ptr<UnivariatePolynomial<Coeff>>,
                const SignCondition&          
        );
        
        /*
         * Constructs a univariate thom encoding with the given polynomial and the sign conditions
         * The polynomial is asserted to be univariate
         */
        ThomEncoding(
                std::shared_ptr<MultivariatePolynomial<Coeff>>,
                const SignCondition&          
        );
                
        /*
         * Constructs a (multivariate) thom encoding with the given polynomial
         * with the given main variable, the sign conditions and the preceding point.
         */
        ThomEncoding(
                std::shared_ptr<MultivariatePolynomial<Coeff>>,
                Variable::Arg, 
                const SignCondition&,               
                std::shared_ptr<ThomEncoding<Coeff>>
        );
        
        /*
         * Constructs a trivial univariate thom encoding for the given rational number.
         */
        ThomEncoding(
                const Coeff& rational,
                Variable::Arg var
        );
        
        /*
         * Constructs a multivariate thom encoding that represents the given rational
         */
        ThomEncoding(
                const Coeff& rational,
                Variable::Arg var,
                std::shared_ptr<ThomEncoding<Coeff>> point
        );
        
        /////////////////////////////////////
        // general public member functions //
        /////////////////////////////////////
        
        /*
         * returns the polynomial
         */
        inline MultivariatePolynomial<Coeff> polynomial() const { return *p; }
        
        /*
         * returns a list of all the underlying polynomials, including the lower levels
         * this is done recursively
         */
        std::vector<MultivariatePolynomial<Coeff>> accumulatePolynomials() const;
        
        /*
         * returns the main variable
         */
        inline Variable::Arg getMainVar() const { return mainVar; }
        
        /*
         * returns the sign condition
         */
        inline SignCondition getSigns() const { return signs; }
        
        /*
         * recursive version of getSigns
         */
        SignCondition accumulateSigns() const;
        
        /*
         * Returns the sign condition of the represented number realized on der(p) = p,p',...,p^{deg(p)}
         */
        SignCondition fullSignCondition() const;
        
        /*
         * Returns the sign of the n-th derivative P^(n) (counting starts from 1)
         */
        Sign operator[](const uint n) const;
        
        /*
         * constructs the (reduced) list of derivatives the internal sign cond. is related to
         * this is precisely the list p',...,p^{deg(p) - 1}
         */
        std::vector<MultivariatePolynomial<Coeff>> reducedDer() const;
        
        /*
         * recursive version of reducedDer()
         */
        std::vector<MultivariatePolynomial<Coeff>> accumulateDer() const;
        
        /*
         * query if the encoding is one or multidimensional
         */
        inline bool isOneDimensional() const { 
                bool res = point.get() == nullptr;
                if(res) assert(polynomial().gatherVariables().size() == 1);
                return res;
        }       
        inline bool isMultiDimensional() const { return !isOneDimensional(); }
        
        /*
         * returns the dimension of the represented point
         * this is done recursively
         */
        uint dimension() const;
                   
        /*
         * checks some basic invariants of a thom encoding object
         */
        bool isConsistent() const;
        
        /*
         * Returns true if and only if the given rational number is encoded by this Thom encoding.
         */
        inline bool represents(const Coeff& rational) const { return *(this) == rational; }
        
        /*
         * returns the sign of the respresented RAN
         */
        Sign sgnReprNum() const;
        
        /*
         * states if the given encoding can be compared by our operators
         */
        template<typename C>
        friend bool areComparable(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs);
        
        ////////////////
        // Comparison //
        ////////////////
        
        /*
         * Compares the RANs represented by the given Thom ecodings
         * This is quite expensive if the underlying polynomials are different.
         * Note that no encodings with different previous points can be compared
         * 
         */
        template<typename C>
        friend ThomComparisonResult compareRational(const ThomEncoding<C>& t, const C& rat);
        template<typename C>
        friend ThomComparisonResult compareThom(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs);
        template<typename C>
        friend ThomComparisonResult compareUnivariate(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs);
        template<typename C>
        friend ThomComparisonResult compareMultivariate(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs);
        
        // compare thom encodings with thom encodings
        template<typename C>
        friend bool operator<(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs) { return compareThom(lhs, rhs) == LESS; }
        template<typename C>
        friend bool operator<=(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs) { return compareThom(lhs, rhs) != GREATER; } 
        template<typename C>
        friend bool operator>(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs) { return compareThom(lhs, rhs) == GREATER; }    
        template<typename C>
        friend bool operator>=(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs) { return compareThom(lhs, rhs) != LESS; } 
        template<typename C>
        friend bool operator==(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs) { return compareThom(lhs, rhs) == EQUAL; }
        template<typename C>
        friend bool operator!=(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs) { return compareThom(lhs, rhs) != EQUAL; }
        
        // comparing with rational numbers
        template<typename C>
        friend bool operator<(const ThomEncoding<C>& lhs, const C& rhs) { return compareRational(lhs, rhs) == LESS; }   
        template<typename C>
        friend bool operator<=(const ThomEncoding<C>& lhs, const C& rhs) { return compareRational(lhs, rhs) != GREATER; }
        template<typename C>
        friend bool operator>(const ThomEncoding<C>& lhs, const C& rhs) { return compareRational(lhs, rhs) == GREATER; }
        template<typename C>
        friend bool operator>=(const ThomEncoding<C>& lhs, const C& rhs) { return compareRational(lhs, rhs) != LESS; }    
        template<typename C>
        friend bool operator==(const ThomEncoding<C>& lhs, const C& rhs) { return compareRational(lhs, rhs) == EQUAL; }  
        template<typename C>
        friend bool operator!=(const ThomEncoding<C>& lhs, const C& rhs) { return compareRational(lhs, rhs) != EQUAL; }
              
        template<typename C>
        friend bool operator<(const C& lhs, const ThomEncoding<C>& rhs) { return rhs > lhs; }      
        template<typename C>
        friend bool operator<=(const C& lhs, const ThomEncoding<C>& rhs) { return rhs >= lhs; }    
        template<typename C>
        friend bool operator>(const C& lhs, const ThomEncoding<C>& rhs) { return rhs < lhs; }   
        template<typename C>
        friend bool operator>=(const C& lhs, const ThomEncoding<C>& rhs) { return rhs <= lhs; }
        template<typename C>
        friend bool operator==(const C& lhs, const ThomEncoding<C>& rhs) { return rhs == lhs; } 
        template<typename C>
        friend bool operator!=(const C& lhs, const ThomEncoding<C>& rhs) { return rhs != lhs; }
             
        /////////////////////////
        // Intermediate points //
        /////////////////////////
        
        // returns a thom encoding representing a number in the interval (lhs, rhs)
        template<typename C, typename Settings>
        friend ThomEncoding<C> intermediatePoint(const ThomEncoding<C>& lhs, const ThomEncoding<C>& rhs);
        
        // adding rational numbers to ran's represented by thom encodings
        template<typename C>
        friend ThomEncoding<C> operator+(const ThomEncoding<C>& rhs, const C& lhs);
        
        template<typename C>
        friend ThomEncoding<C> operator+(const C& lhs, const ThomEncoding<C>& rhs) { return rhs + lhs; }
            
        /*
         * output operator
         */
        template <typename C, typename Settings = ThomDefaultSettings>
	friend std::ostream& operator<<(std::ostream& os, const ThomEncoding<C>& rhs) {
                os << "(" << rhs.polynomial() << " with mainVar " << rhs.getMainVar() << ", " << rhs.fullSignCondition() << ")" << std::endl;
                if(rhs.isMultiDimensional()) {
                        os << *(rhs.point);
                }
                return os;
        }
};


#include "ThomEncoding.tpp"
} // namespace carl



