#pragma once

#include "../config.h"
#include "PrimitiveEuclidean.h"
#include "../MultivariatePolynomial.h"
#include "../../numbers/typetraits.h"

#include "../../converter/CoCoAAdaptor.h"
#include "../../converter/OldGinacConverter.h"

namespace carl {

namespace gcd_detail {
	template<typename Polynomial>
	Variable select_variable(const Polynomial& p1, const Polynomial& p2) {
		auto v1 = carl::variables(p1).as_vector();
		auto v2 = carl::variables(p2).as_vector();
		std::vector<Variable> common;
		
		std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(), std::inserter(common,common.begin()));
		if (common.empty()) {
			return Variable::NO_VARIABLE;
		} else {
			return *common.begin();
		}
	}

	template<typename Polynomial>
	Polynomial gcd_calculate(const Polynomial& a, const Polynomial& b) {
		Variable x = select_variable(a, b);
		if (x == Variable::NO_VARIABLE) {
			return Polynomial(1);
		}
		CARL_LOG_INEFFICIENT();
		auto A = a.toUnivariatePolynomial(x);
		auto B = b.toUnivariatePolynomial(x);
		Polynomial result(primitive_euclidean(A.normalized(), B.normalized()));
		if (carl::isNegative(result.lcoeff()) && !(carl::isNegative(a.lcoeff()) && carl::isNegative(b.lcoeff()))) {
			return -result;
		}
		return result;
	}
}

template<typename C, typename O, typename P>
MultivariatePolynomial<C,O,P> gcd(const MultivariatePolynomial<C,O,P>& a, const MultivariatePolynomial<C,O,P>& b) {
	CARL_LOG_DEBUG("carl.core.gcd", "gcd(" << a << ", " << b << ")");
	assert(!isZero(a));
	assert(!isZero(b));

	if (isOne(a) || isOne(b)) {
		return MultivariatePolynomial<C,O,P>(1);
	}
	if (a.isConstant() && b.isConstant()) {
		CARL_LOG_DEBUG("carl.core.gcd", "gcd(" << a << ", " << b << ") = " << carl::gcd(a.constantPart(), b.constantPart()));
		return MultivariatePolynomial<C,O,P>(carl::gcd(a.constantPart(), b.constantPart()));
	}
	if (a.isConstant() || b.isConstant()) {
		return MultivariatePolynomial<C,O,P>(1);
	}

	auto s = overloaded {
	#if defined USE_GINAC
		[](const MultivariatePolynomial<cln::cl_RA,O,P>& n1, const MultivariatePolynomial<cln::cl_RA,O,P>& n2){ return ginacGcd<MultivariatePolynomial<cln::cl_RA,O,P>>( n1, n2 ); },
		[](const MultivariatePolynomial<cln::cl_I,O,P>& n1, const MultivariatePolynomial<cln::cl_I,O,P>& n2){ return ginacGcd<MultivariatePolynomial<cln::cl_I,O,P>>( n1, n2 ); },
	#endif
	#if defined USE_COCOA
		[](const MultivariatePolynomial<mpq_class,O,P>& n1, const MultivariatePolynomial<mpq_class,O,P>& n2){ CoCoAAdaptor<MultivariatePolynomial<mpq_class,O,P>> c({n1, n2}); return c.gcd(n1,n2); },
		[](const MultivariatePolynomial<mpz_class,O,P>& n1, const MultivariatePolynomial<mpz_class,O,P>& n2){ CoCoAAdaptor<MultivariatePolynomial<mpz_class,O,P>> c({n1, n2}); return c.gcd(n1,n2); }
	#else
		[](const MultivariatePolynomial<mpq_class,O,P>& n1, const MultivariatePolynomial<mpq_class,O,P>& n2){ return gcd_detail::gcd_calculate(n1,n2); },
		[](const MultivariatePolynomial<mpz_class,O,P>& n1, const MultivariatePolynomial<mpz_class,O,P>& n2){ return gcd_detail::gcd_calculate(n1,n2); }
	#endif
	};
	CARL_LOG_DEBUG("carl.core.gcd", "gcd(" << a << ", " << b << ")");
	auto res = s(a, b);
	CARL_LOG_DEBUG("carl.core.gcd", "gcd(" << a << ", " << b << ") = " << res);
	return res;
}
}