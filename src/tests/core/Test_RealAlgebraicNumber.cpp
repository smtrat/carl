#include "gtest/gtest.h"
#include <map>

#include "carl/core/UnivariatePolynomial.h"
#include "carl/ran/ran.h"

#include "../Common.h"

using namespace carl;


TEST(RealAlgebraicNumber, Evaluation)
{
    Variable y = freshRealVariable("skoY");
    Variable x = freshRealVariable("skoX");

	MultivariatePolynomial<Rational> mpx(x);
	MultivariatePolynomial<Rational> mpy(y);

	MultivariatePolynomial<Rational> lin(Rational(-72)*carl::pow(mpx, 6) + Rational(3024)*carl::pow(mpx, 4) + Rational(-60480)*carl::pow(mpx, 2));
	UnivariatePolynomial<MultivariatePolynomial<Rational>> p(y, {
		MultivariatePolynomial<Rational>(0),
		lin,
		MultivariatePolynomial<Rational>(0),
		MultivariatePolynomial<Rational>(60480),
		MultivariatePolynomial<Rational>(0),
		MultivariatePolynomial<Rational>(-3024),
		MultivariatePolynomial<Rational>(0),
		MultivariatePolynomial<Rational>(72),
		MultivariatePolynomial<Rational>(0),
		MultivariatePolynomial<Rational>(-1)
	});

	//UnivariatePolynomial<Rational> px(y, {-3528000, 0, 882000, 0, -74760, 0, 3444, 0, -84, 0, 1});
	UnivariatePolynomial<Rational> px(y, std::initializer_list<Rational>{-3528000, 0, 882000, 0, -74760, 0, 3444, 0, -84, 0, 1});
	Interval<Rational> ix(Rational("2927288666429")/Rational("1099511627776"), BoundType::STRICT, Rational("1463644333215")/Rational("549755813888"), BoundType::STRICT);
	UnivariatePolynomial<Rational> py(y, std::initializer_list<Rational>{Rational("-18289152000"), 0, Rational("4572288000"), 0, Rational("-387555840"), 0, 18156096, 0, -556416, 0, 11232, 0, -144, 0, 1});
	Interval<Rational> iy(Rational(-147580509)/24822226, BoundType::STRICT, Rational(-73113831)/12411113, BoundType::STRICT);

	/* m = {
	 * skoY : (IR ]-147580509/24822226, -73113831/12411113[, (1)*skoY^14 + (-144)*skoY^12 + (11232)*skoY^10 + (-556416)*skoY^8 + (18156096)*skoY^6 + (-387555840)*skoY^4 + (4572288000)*skoY^2 + -18289152000),
	 * skoX : (IR ]85/32, 341/128[, (1)*skoX^10 + (-84)*skoX^8 + (3444)*skoX^6 + (-74760)*skoX^4 + (882000)*skoX^2 + -3528000)
	 * }
	 */
}



TEST(RealAlgebraicNumber, EvalBug)
{
	// Evaluating skoY^2 + skoX^2 <= 0 on {skoY : (IR ]-212079/131072, -1696631/1048576[, __r^2 + __r^1 + -1), skoD : (NR 1), skoX : (NR 0)} -> 1

    Variable y = freshRealVariable("skoY");
    Variable x = freshRealVariable("skoX");
	MultivariatePolynomial<Rational> mpx(x);
	MultivariatePolynomial<Rational> mpy(y);
	MultivariatePolynomial<Rational> poly(carl::pow(mpx, 2) + carl::pow(mpy, 2));
	Constraint<MultivariatePolynomial<Rational>> constr(poly,carl::Relation::LEQ);

	Variable h = freshRealVariable("h"); 
	UnivariatePolynomial<Rational> py(h, std::initializer_list<Rational>{-1, 1, 1});
	Interval<Rational> iy(Rational("-212079")/Rational("131072"), BoundType::STRICT, -Rational("1696631")/Rational("1048576"), BoundType::STRICT);
	RealAlgebraicNumber<Rational> ry = RealAlgebraicNumber<Rational>::create_safe(py, iy);

	carl::ran::RANMap<Rational> eval;
	eval.emplace(y,ry);
	eval.emplace(x,Rational(0));

	auto res = carl::evaluate(constr, eval);
	EXPECT_FALSE((bool)res);
}

TEST(RealAlgebraicNumber, EvalBug2)
{
	// Poly : 2*x1 + 3*x2  
	// Assignment : {x1 : (NR -4), x2 : (IR ]-3, -2[, (2)*__r^3 + (6)*__r^2 + (4)*__r^1 + 11), x0 : (NR -8)} 
	Variable x0 = freshRealVariable("x0");
	Variable x1 = freshRealVariable("x1");
	Variable x2 = freshRealVariable("x2");
	MultivariatePolynomial<Rational> mpx1(x1) ;
	MultivariatePolynomial<Rational> mpx2(x2) ;

	//build poly
	MultivariatePolynomial<Rational> poly(Rational("2")*mpx1 + Rational("3")*mpx2) ;

	//build ran
	Variable h = freshRealVariable("h"); 
	UnivariatePolynomial<Rational> py(h, std::initializer_list<Rational>{11, 4, 6, 2});
	Interval<Rational> iy(Rational("-3"), BoundType::STRICT, Rational("-2"), BoundType::STRICT);
	RealAlgebraicNumber<Rational> ry = RealAlgebraicNumber<Rational>::create_safe(py, iy);

	//build assignment
	carl::ran::RANMap<Rational> eval;
	eval.emplace(x0, Rational("-8"));
	eval.emplace(x1, Rational("-4"));
	eval.emplace(x2, ry);

	auto res = carl::evaluate(poly, eval);
	EXPECT_TRUE((bool) res);
}




