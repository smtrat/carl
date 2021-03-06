#pragma once

#include "ThomEncoding.h"
#include "../ran_operations.h"


#include <memory>

namespace carl {

// TODO adapt to new interface

template<typename Number>
struct real_algebraic_number_thom {
	template<typename Num>
	friend bool operator==(const real_algebraic_number_thom<Num>& lhs, const real_algebraic_number_thom<Num>& rhs);
	template<typename Num>
	friend bool operator<(const real_algebraic_number_thom<Num>& lhs, const real_algebraic_number_thom<Num>& rhs);
private:
	struct Content {
		ThomEncoding<Number> te;

		Content(const ThomEncoding<Number>& t):
			te(t)
		{}
	};
	std::shared_ptr<Content> mContent;
public:
	real_algebraic_number_thom(const ThomEncoding<Number>& te):
		mContent(std::make_shared<Content>(te))
	{}

	auto& thom_encoding() {
		return mContent->te;
	}
	const auto& thom_encoding() const {
		return mContent->te;
	}

	const auto& polynomial() const {
		return thom_encoding().polynomial();
	}
	const auto& mainVar() const {
		return thom_encoding().mainVar();
	}
	auto sign_condition() const {
		return thom_encoding().signCondition();
	}
	const auto& point() const {
		return thom_encoding().point();
	}

	std::size_t size() const {
		return thom_encoding().dimension();
	}

	std::size_t dimension() const {
		return thom_encoding().dimension();
	}

	bool is_integral() const {
		return thom_encoding().is_integral();
	}
	bool is_zero() const {
		return thom_encoding().is_zero();
	}
	bool contained_in(const Interval<Number>& i) const {
		return thom_encoding().containedIn(i);
	}

	Number integer_below() const {
		return thom_encoding().integer_below();
	}
	Sign sgn() const {
		return thom_encoding().sgn();
	}
	Sign sgn(const UnivariatePolynomial<Number>& p) const {
		return thom_encoding().sgn(p);
	}
};

template<typename Number>
Number branching_point(const real_algebraic_number_thom<Number>& n) {
	return n.thom_encoding().get_number();
}

template<typename Number>
Number evaluate(const MultivariatePolynomial<Number>& p, std::map<Variable, real_algebraic_number_thom<Number>>& m) {
	//using Polynomial = MultivariatePolynomial<Number>;
	
	CARL_LOG_INFO("carl.ran.thom",
			"\n****************************\n"
			<< "Thom evaluate\n"
			<< "****************************\n"
			<< "p = " << p << "\n"
			<< "m = " << m << "\n"
			<< "****************************\n");
	for(const auto& entry : m) { 
		assert(entry.first == entry.second.thom_encoding().mainVar());
	}
	assert(m.size() > 0);
	
	std::map<Variable, real_algebraic_number_thom<Number>>& m_prime(m);
	auto it = m_prime.begin();
	while(it != m_prime.end()) {
			if(!p.has(it->first)) {
					CARL_LOG_TRACE("carl.thom.evaluation", "removing " << it->first);
					it = m_prime.erase(it);
			} else {
					it++;
			}
	}
	
	std::map<Variable, ThomEncoding<Number>> mTE;
	for(const auto& entry : m_prime) {
			mTE.insert(std::make_pair(entry.first, entry.second.thom_encoding()));
	}
	
	CARL_LOG_ASSERT("carl.thom.evaluation", p.gatherVariables().size() == mTE.size(), "p = " << p << ", mTE = " << mTE);
	
	if(mTE.size() == 1) {
			int sgn = int(mTE.begin()->second.signOnPolynomial(p));
			CARL_LOG_TRACE("carl.thom.evaluation", "sign of evaluated polynomial is " << sgn);
			return Number(sgn);
	}
	
	CARL_LOG_TRACE("carl.thom.evaluation", "mTE = " << mTE);
	
	ThomEncoding<Number> point = ThomEncoding<Number>::analyzeTEMap(mTE);
	int sgn = int(point.signOnPolynomial(p));
	CARL_LOG_TRACE("carl.thom.", "sign of evaluated polynomial is " << sgn);
	return Number(sgn);
	
}

template<typename Number, typename Poly>
bool evaluate(const Constraint<Poly>& c, std::map<Variable, real_algebraic_number_thom<Number>>& m) {
	auto res = evaluate(c.lhs(), m);
	return evaluate(res, c.relation());
}

template<typename Number>
real_algebraic_number_thom<Number> abs(const real_algebraic_number_thom<Number>& n) {
	assert(false);
	return n;
}

template<typename Number>
real_algebraic_number_thom<Number> sample_above(const real_algebraic_number_thom<Number>& n) {
	return n.thom_encoding() + Number(1);
}
template<typename Number>
real_algebraic_number_thom<Number> sample_below(const real_algebraic_number_thom<Number>& n) {
	return n.thom_encoding() + Number(-1);
}
template<typename Number>
real_algebraic_number_thom<Number> sample_between(const real_algebraic_number_thom<Number>& lower, const real_algebraic_number_thom<Number>& upper) {
	return ThomEncoding<Number>::intermediatePoint(lower.thom_encoding(), upper.thom_encoding());
}
template<typename Number>
Number sample_between(const real_algebraic_number_thom<Number>& lower, const Number& upper) {
	return ThomEncoding<Number>::intermediatePoint(lower.thom_encoding(), upper);
}
template<typename Number>
Number sample_between(const Number& lower, const real_algebraic_number_thom<Number>& upper) {
	return ThomEncoding<Number>::intermediatePoint(lower, upper.thom_encoding());
}

template<typename Number>
Number floor(const real_algebraic_number_thom<Number>& n) {
	return carl::floor(get_interval(n).lower());
}
template<typename Number>
Number ceil(const real_algebraic_number_thom<Number>& n) {
	return carl::ceil(get_interval(n).upper());
}

template<typename Number>
bool operator==(const real_algebraic_number_thom<Number>& lhs, const real_algebraic_number_thom<Number>& rhs) {
	if (lhs.mContent.get() == rhs.mContent.get()) return true;
	return lhs.thom_encoding() == rhs.thom_encoding();
}

template<typename Number>
bool operator==(const real_algebraic_number_thom<Number>& lhs, const Number& rhs) {
	return lhs.thom_encoding() == rhs;
}

template<typename Number>
bool operator==(const Number& lhs, const real_algebraic_number_thom<Number>& rhs) {
	return lhs == rhs.thom_encoding();
}

template<typename Number>
bool operator<(const real_algebraic_number_thom<Number>& lhs, const real_algebraic_number_thom<Number>& rhs) {
	if (lhs.mContent.get() == rhs.mContent.get()) return false;
	return lhs.thom_encoding() < rhs.thom_encoding();
}

template<typename Number>
bool operator<(const real_algebraic_number_thom<Number>& lhs, const Number& rhs) {
	return lhs.thom_encoding() == rhs;
}

template<typename Number>
bool operator<(const Number& lhs, const real_algebraic_number_thom<Number>& rhs) {
	return lhs == rhs.thom_encoding();
}

template<typename Num>
std::ostream& operator<<(std::ostream& os, const real_algebraic_number_thom<Num>& rhs) {
	os << "(TE " << rhs.polynomial() << " in " << rhs.mainVar() << ", " << rhs.sign_condition();
	if (rhs.dimension() > 1) {
		os << " OVER " << rhs.point();
	}
	os << ")";
	return os;
}

template<typename Number>
struct is_ran<real_algebraic_number_thom<Number>> { 
  static const bool value = true;
};

}


namespace std {
template<typename Number>
struct hash<carl::real_algebraic_number_z3<Number>> {
    std::size_t operator()(const carl::real_algebraic_number_z3<Number>& n) const {
		return carl::hash_all(n.integer_below());
	}
};
}