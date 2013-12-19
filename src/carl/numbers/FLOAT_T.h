/**
 * General class for floating point numbers with different formats. Extend to
 * other types if necessary.
 * 
 * @file FLOAT_T.h
 * @author  Stefan Schupp
 * @since   2013-10-14
 * @version 2013-10-14
 */

#pragma once

#include <string>
#include <iostream>
#include <assert.h>

#include "config.h"
#include <math.h>
#ifdef USE_MPFR_FLOAT
#include <mpfr.h>
#endif

namespace carl
{
		// nearest, towards zero, towards infty, towards -infty, away from zero
//    enum CARL_RND
//    {
//        CARL_RNDN=0,
//        CARL_RNDZ=1,
//        CARL_RNDU=2,
//        CARL_RNDD=3,
//        CARL_RNDA=4
//    };

    typedef unsigned long precision;
    
    template<typename FloatType>
    class FLOAT_T;
    
    template<typename T1, typename T2>
    struct FloatConv
    {
        FLOAT_T<T1> operator() ( const FLOAT_T<T2>& _op2 ) const
        {
            return FLOAT_T<T1>(_op2.toDouble());
        }
    };
    
    template<typename FloatType>
    class FLOAT_T
    {
        private:
            FloatType mValue;
        
        public:
    	
            /**
             * Constructors & Destructors
             */

            FLOAT_T<FloatType>()
            {
				assert(is_float(FloatType));
                mValue = 0;
            }
        
            FLOAT_T<FloatType>(const double _double)
            {
				assert(is_float(FloatType));
                mValue = _double;
            }
            
            FLOAT_T<FloatType>(const float _float)
            {
				assert(is_float(FloatType));
                mValue = _float;
            }
            
            FLOAT_T<FloatType>(const int _int)
            {
				assert(is_float(FloatType));
                mValue = _int;
            }
            
            FLOAT_T<FloatType>(const FLOAT_T<FloatType>& _float) : mValue(_float.mValue)
            {
				assert(is_float(FloatType));
			}
            
            ~FLOAT_T()
            {}

            
            /*******************
            * Getter & Setter *
            *******************/

           const FloatType& getValue() const
           {
               return mValue;
           }

           precision getPrec() const
           {
               // TODO
               return 0;
           }

           FLOAT_T<FloatType>& setPrec( const precision& _prec)
           {
               // TODO
               return *this;
           }
            
            /*************
             * Operators *
             *************/
            
            FLOAT_T<FloatType>& operator = (const FLOAT_T<FloatType> & _rhs)
            {
                mValue = _rhs.mValue;
                return *this;
            }
            
            /**
             * Boolean operators 
             */
            
            bool operator == ( const FLOAT_T<FloatType>& _rhs) const
            {
                return mValue == _rhs.mValue;
            }

            bool operator != ( const FLOAT_T<FloatType> & _rhs) const
            {
                return mValue != _rhs.mValue;
            }

            bool operator > ( const FLOAT_T<FloatType> & _rhs) const
            {
                return mValue > _rhs.mValue;
            }

            bool operator < ( const FLOAT_T<FloatType> & _rhs) const
            {
                return mValue < _rhs.mValue;
            }

            bool operator <= ( const FLOAT_T<FloatType> & _rhs) const
            {
                return mValue <= _rhs.mValue;
            }

            bool operator >= ( const FLOAT_T<FloatType> & _rhs) const
            {
                return mValue >= _rhs.mValue;
            }
            
            /**
             * arithmetic operations
             */
            
            FLOAT_T& add_assign( const FLOAT_T<FloatType>& _op2, CARL_RND _rnd = CARL_RNDN)
            {
                // TODO: Include rounding
                mValue = mValue + _op2.mValue;
                return *this;
            }
            
            void add( FLOAT_T& _result, const FLOAT_T& _op2, CARL_RND _rnd = CARL_RNDN) const
            {
                // TODO: Include rounding
                _result.mValue = mValue + _op2.mValue;
            }

            FLOAT_T& sub_assign(const FLOAT_T& _op2, CARL_RND _rnd = CARL_RNDN)
            {
                // TODO: Include rounding
                mValue = mValue - _op2.mValue;
                return *this;
            }
            
            void sub(FLOAT_T& _result, const FLOAT_T& _op2, CARL_RND _rnd = CARL_RNDN) const
            {
                // TODO: Include rounding
                _result.mValue = mValue - _op2.mValue;
            }
            
            FLOAT_T& mul_assign(const FLOAT_T& _op2, CARL_RND _rnd = CARL_RNDN)
            {
                // TODO: Include rounding
                mValue = mValue * _op2.mValue;
                // mValue.mul_assign(_op2.mValue, _rnd);
                return *this;
            }
            
            void mul(FLOAT_T& _result, const FLOAT_T& _op2, CARL_RND _rnd = CARL_RNDN) const
            {
                // TODO: Include rounding
                _result.mValue = mValue * _op2.mValue;
            }
            
            FLOAT_T& div_assign(const FLOAT_T& _op2, CARL_RND _rnd = CARL_RNDN)
            {
                assert( _op2 != 0 );
                // TODO: Include rounding
                mValue = mValue / _op2.mValue;
                return *this;
            }
            
            void div(FLOAT_T& _result, const FLOAT_T& _op2, CARL_RND _rnd = CARL_RNDN) const
            {
                assert( _op2 != 0 );
                // TODO: Include rounding
                _result.mValue = mValue / _op2.mValue;
            }
            
            /**
             * special operators
             */

            FLOAT_T& sqrt_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = sqrt(mValue);
                return *this;
            }
            
            void sqrt(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
            	_result.mValue = sqrt(mValue);
            }

            FLOAT_T& cbrt_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = cbrt(mValue);
                return *this;
            }
            
            FLOAT_T& cbrt(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
            	_result.mValue = cbrt(mValue);
            }

            FLOAT_T& root_assign(unsigned long int _k, CARL_RND _rnd = CARL_RNDN)
            {
                // TODO
                return *this;
            }
            
            void root(FLOAT_T& _result, unsigned long int _k, CARL_RND _rnd = CARL_RNDN)
            {
            	// TODO
            }

            FLOAT_T& pow_assign(unsigned long int _exp, CARL_RND _rnd = CARL_RNDN)
            {
                mValue = pow(mValue, _exp);
                return *this;
            }
            
            void pow(FLOAT_T& _result, unsigned long int _exp, CARL_RND _rnd = CARL_RNDN)
            {
            	_result.mValue = pow(mValue, _exp);
            }

            FLOAT_T& abs_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = abs(mValue);
                return *this;
            }
            
            void abs(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
            	_result.mValue = abs(mValue);
            }
            
            FLOAT_T& exp_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = exp(mValue);
                return *this;
            }
            
            void exp(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = exp(mValue);
            }
            
            FLOAT_T& sin_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = sin(mValue);
                return *this;
            }
            
            void sin(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = sin(mValue);
            }
            
            FLOAT_T& cos_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = cos(mValue);
                return *this;
            }
            
            void cos(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = cos(mValue);
            }
            
            FLOAT_T& log_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = log(mValue);
                return *this;
            }
            
            void log(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = log(mValue);
            }
            
            FLOAT_T& tan_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = tan(mValue);
                return *this;
            }
            
            void tan(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = tan(mValue);
            }
            
            FLOAT_T& asin_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = asin(mValue);
                return *this;
            }
            
            void asin(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = asin(mValue);
            }
            
            FLOAT_T& acos_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = acos(mValue);
                return *this;
            }
            
            void acos(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = acos(mValue);
            }
            
            FLOAT_T& atan_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = atan(mValue);
                return *this;
            }
            
            void atan(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = atan(mValue);
            }

            FLOAT_T& sinh_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = sinh(mValue);
                return *this;
            }
            
            void sinh(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = sinh(mValue);
            }
            
            FLOAT_T& cosh_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = cosh(mValue);
                return *this;
            }
            
            void cosh(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = cosh(mValue);
            }
            
            FLOAT_T& tanh_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = tanh(mValue);
                return *this;
            }
            
            void tanh(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = tanh(mValue);
            }
            
            FLOAT_T& asinh_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = asinh(mValue);
                return *this;
            }
            
            void asinh(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = asinh(mValue);
            }
            
            FLOAT_T& acosh_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = acosh(mValue);
                return *this;
            }
            
            void acosh(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = acosh(mValue);
            }
            
            FLOAT_T& atanh_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = atanh(mValue);
                return *this;
            }
            
            void atanh(FLOAT_T& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result.mValue = atanh(mValue);
            }
            
            void floor(int& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result = floor(mValue);
            }
            
            FLOAT_T& floor_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = floor(mValue);
                return *this;
            }
            
            void ceil(int& _result, CARL_RND _rnd = CARL_RNDN)
            {
                _result = ceil(mValue);
            }
            
            FLOAT_T& ceil_assign(CARL_RND _rnd = CARL_RNDN)
            {
                mValue = ceil(mValue);
                return *this;
            }
            
            /**
             * conversion operators
             */
            double toDouble(CARL_RND _rnd=CARL_RND::CARL_RNDN) const
            {
                return (double) mValue;
            }
            

            friend std::ostream & operator<< (std::ostream& ostr, const FLOAT_T<FloatType>& p) {
                ostr << p.toString();
                return ostr;
            }
            
            friend bool operator== (const FLOAT_T<FloatType>& _lhs, const int _rhs)
            {
                return _lhs.mValue == _rhs;
            }
            
            friend bool operator== (const int _lhs, const FLOAT_T<FloatType>& _rhs)
            {
                return _rhs == _lhs;
            }
            
            friend bool operator== (const FLOAT_T<FloatType>& _lhs, const double _rhs)
            {
                return _lhs.mValue == _rhs;
            }
            
            friend bool operator== (const double _lhs, const FLOAT_T<FloatType>& _rhs)
            {
                return _rhs == _lhs;
            }
            
            friend bool operator== (const FLOAT_T<FloatType>& _lhs, const float _rhs)
            {
                return _lhs.mValue == _rhs;
            }
            
            friend bool operator== (const float _lhs, const FLOAT_T<FloatType>& _rhs)
            {
                return _rhs == _lhs;
            }
            
            /**
             * Operators
             */
            
            friend FLOAT_T<FloatType> operator +(const FLOAT_T<FloatType>& _lhs, const FLOAT_T<FloatType>& _rhs)
            {
                return FLOAT_T<FloatType>(_lhs.mValue + _rhs.mValue);
            }
		
			friend FLOAT_T<FloatType> operator +(const FLOAT_T<FloatType>& _lhs, const FloatType& _rhs)
			{
				return FLOAT_T<FloatType>(_lhs.mValue + _rhs);
			}
        
			friend FLOAT_T<FloatType> operator +(const FloatType& _lhs, const FLOAT_T<FloatType>& _rhs)
			{
				return _rhs + _lhs;
			}
		
            friend FLOAT_T<FloatType> operator -(const FLOAT_T<FloatType>& _lhs, const FLOAT_T<FloatType>& _rhs)
            {
                return FLOAT_T<FloatType>(_lhs.mValue - _rhs.mValue);
            }

			friend FLOAT_T<FloatType> operator -(const FLOAT_T<FloatType>& _lhs, const FloatType& _rhs)
			{
				return FLOAT_T<FloatType>(_lhs.mValue - _rhs);
			}
			
			friend FLOAT_T<FloatType> operator -(const FloatType& _lhs, const FLOAT_T<FloatType>& _rhs)
			{
				return FLOAT_T<FloatType>(_lhs - _rhs.mValue);
			}
		
            friend FLOAT_T<FloatType> operator *(const FLOAT_T<FloatType>& _lhs, const FLOAT_T<FloatType>& _rhs)
            {
                return FLOAT_T<FloatType>(_lhs.mValue * _rhs.mValue);
            }
			
			friend FLOAT_T<FloatType> operator *(const FloatType& _lhs, const FLOAT_T<FloatType>& _rhs)
			{
				return FLOAT_T<FloatType>(_lhs * _rhs.mValue);
			}
			
			friend FLOAT_T<FloatType> operator *(const FLOAT_T<FloatType>& _lhs, const FloatType& _rhs)
			{
				return FLOAT_T<FloatType>(_lhs.mValue * _rhs);
			}
        
            friend FLOAT_T<FloatType> operator /(const FLOAT_T<FloatType>& _lhs, const FLOAT_T<FloatType>& _rhs)
            {
                return FLOAT_T<FloatType>(_lhs.mValue / _rhs.mValue);
            }

			friend FLOAT_T<FloatType> operator /(const FLOAT_T<FloatType>& _lhs, const FloatType& _rhs)
			{
				return FLOAT_T<FloatType>(_lhs.mValue / _rhs);
			}
			
			friend FLOAT_T<FloatType> operator /(const FloatType& _lhs, const FLOAT_T<FloatType>& _rhs)
			{
				return FLOAT_T<FloatType>(_lhs / _rhs.mValue);
			}
            friend FLOAT_T<FloatType>& operator ++(FLOAT_T<FloatType>& _num)
            {
                _num.mValue +=1;
                return _num;
            }
            
            friend FLOAT_T<FloatType>& operator --(FLOAT_T<FloatType>& _num)
            {
                _num.mValue -=1;
                return _num;
            }
            
            FLOAT_T<FloatType>& operator +=(const FLOAT_T<FloatType>& _rhs)
            {
                mValue + _rhs.mValue;
                return *this;
            }
        
			FLOAT_T<FloatType>& operator +=(const FloatType& _rhs)
			{
				mValue + _rhs;
				return *this;
			}
		
            FLOAT_T<FloatType>& operator -=(const FLOAT_T<FloatType>& _rhs)
            {
                mValue - _rhs.mValue;
                return *this;
            }
        
			FLOAT_T<FloatType>& operator -=(const FloatType& _rhs)
			{
				mValue - _rhs;
				return *this;
			}
		
            FLOAT_T<FloatType>& operator *=(const FLOAT_T<FloatType>& _rhs)
            {
                mValue * _rhs.mValue;
                return *this;
            }
        
			FLOAT_T<FloatType>& operator *=(const FloatType& _rhs)
			{
				mValue * _rhs;
				return *this;
			}
		
            FLOAT_T<FloatType>& operator /=(const FLOAT_T<FloatType>& _rhs)
            {
                mValue / _rhs.mValue;
                return *this;
            }
        
			FLOAT_T<FloatType>& operator /=(const FloatType& _rhs)
			{
				mValue / _rhs;
				return *this;
			}
		
            /**
             * Auxiliary Functions
             */
            
            std::string toString() const
            {
//                std::stringstream str;
//                str << mValue;   
//                return str.str();
                return std::to_string(mValue);
            }
    };
#ifdef USE_MPFR_FLOAT
    #include "floatTypes/mpfr_float.tpp"
#endif

}
