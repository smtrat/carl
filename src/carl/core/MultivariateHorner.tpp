 /**
 * @file	MultivariateHorner.tpp
 * @author	Lukas Netz
 *
 */

//#define DEBUG_HORNER 

namespace carl
{
	//Constructor for Greedy I
	template< typename PolynomialType, Strategy Strat> 
	MultivariateHorner< PolynomialType, Strat>::MultivariateHorner (const PolynomialType&& inPut) {
	#ifdef DEBUG_HORNER
		std::cout << __func__ << " (GreedyI constr) P: " << inPut << std::endl;
	#endif
	static_assert(!(Strat==GREEDY_II)&&!(Strat==GREEDY_IIs), "Strategy requires Interval map");

	Interval<CoeffType> dummy(0);
	Variable dummyVar = Variable::NO_VARIABLE;
	std::map<Variable, Interval<CoeffType>> map;

	map[dummyVar] = dummy;


	//Create Horner Scheme Recursivly
	MultivariateHorner< PolynomialType, Strat > root (std::move(inPut), Strat, map );
	
 	//Part after recursion
 	if (Strat == GREEDY_Is || Strat == GREEDY_IIs)
 	{
 		auto root_ptr =std::make_shared<MultivariateHorner< PolynomialType, Strat > >(root);
 		root_ptr = simplify( root_ptr );	
 		root = *root_ptr;
 	}

 	//Apply all changes
 	*this = root;
	};

	//Constructor for Greedy II and Greedy I
	template< typename PolynomialType, Strategy Strat> 
	MultivariateHorner< PolynomialType, Strat>::MultivariateHorner (const PolynomialType&& inPut, std::map<Variable, Interval<CoeffType>>& map) {
	#ifdef DEBUG_HORNER
		std::cout << __func__ << " (GreedyII constr) P: " << inPut << std::endl;
	#endif

	//Create Horner Scheme Recursivly
	MultivariateHorner< PolynomialType, Strat > root (std::move(inPut), Strat, map);
	
 	//Part after recursion
 	if (Strat == GREEDY_Is || Strat == GREEDY_IIs)
 	{
 		auto root_ptr =std::make_shared<MultivariateHorner< PolynomialType, Strat > >(root);
 		root_ptr = simplify( root_ptr );	
 		root = *root_ptr;
 	}

 	//Apply all changes
 	*this = root;
	};


	//Constructor for Greedy I/II creates recursive Datastruckture
	template< typename PolynomialType, Strategy Strat >
	MultivariateHorner< PolynomialType, Strat>::MultivariateHorner (const PolynomialType&& inPut, Strategy s, std::map<Variable, Interval<CoeffType>>& map) 
	{


		std::set<Variable>::iterator variableIt;
		std::set<Variable>::iterator selectedVariable;
		std::set<Variable> allVariablesinPolynome;
		inPut.gatherVariables(allVariablesinPolynome);

		Interval<CoeffType> currentInterval(0);
		CoeffType delta = constant_zero<CoeffType>::get();
		CoeffType bestDelta = constant_zero<CoeffType>::get();

		unsigned int monomials_containingChoosenVar = 0;

		if (allVariablesinPolynome.size() != 0)
		{
			//Detecting amounts of Variables in Monomials
			for (variableIt = allVariablesinPolynome.begin(); variableIt != allVariablesinPolynome.end(); variableIt++)
			{	
				if (s == GREEDY_I || s == GREEDY_Is)
				{
					unsigned int monomialCounter = 0;
					typename PolynomialType::TermsType::const_iterator polynomialIt;
					for (polynomialIt = inPut.begin(); polynomialIt != inPut.end(); polynomialIt++)
					{
						if (polynomialIt->has(*variableIt))
						{
							monomialCounter++;
						}
					}
					
					//saving most promising Variable for later 
					if ( monomialCounter >= monomials_containingChoosenVar ) {
						monomials_containingChoosenVar = monomialCounter;
						selectedVariable = variableIt;					
					}
				}

				if (s == GREEDY_II || s == GREEDY_IIs)
				{
					typename PolynomialType::TermsType::const_iterator polynomialIt;
					CoeffType accMonomEval = constant_zero<CoeffType>::get();
					CoeffType accMonomDivEval = constant_zero<CoeffType>::get();

					for (polynomialIt = inPut.begin(); polynomialIt != inPut.end(); polynomialIt++)
					{
						if (polynomialIt->has(*variableIt))
						{						
							Term< typename PolynomialType::CoeffType > currentTerm = *polynomialIt;
							Term< typename PolynomialType::CoeffType > currentTerm_div = *polynomialIt;
							
							currentTerm_div.divide(*variableIt, currentTerm_div);

							currentInterval = IntervalEvaluation::evaluate( currentTerm_div, map );

							accMonomEval += IntervalEvaluation::evaluate( currentTerm, map ).diameter();
							accMonomDivEval += currentInterval.diameter();

							//std::cout << "[!]cT: "<< currentTerm << " | " << IntervalEvaluation::evaluate( currentTerm, map ).diameter() << ", " << accMonomEval << " D: " << currentInterval.diameter() << "," << accMonomDivEval << std::endl;
						}
					}		
					accMonomDivEval *= map.find(*variableIt)->second.diameter();
					delta = accMonomDivEval - accMonomEval;
					
					if (delta > bestDelta )
					{
						#ifdef DEBUG_HORNER
							std::cout << "update Delta...";
						#endif

						bestDelta = delta;
						selectedVariable = variableIt;
					}
					//std::cout << *variableIt << " D: "<< delta  << " BD: "<< bestDelta << "\n" << std::endl;
	 			}			
			}

			//Setting the choosen Variable for the current Hornerscheme iterartion

			if (*selectedVariable == NULL || bestDelta == constant_zero<CoeffType>::get() )
			{
				selectedVariable = allVariablesinPolynome.begin();
			}
			this->setVariable(*selectedVariable); 

			#ifdef DEBUG_HORNER
			std::cout << __func__ << "    Polynome: " << inPut << std::endl;
			//std::cout << "Choosen Var: " << *selectedVariable << std::endl;
			#endif
			
			typename PolynomialType::TermsType::const_iterator polynomialIt;
			typename PolynomialType::TermType tempTerm;

			typename PolynomialType::PolyType h_independentPart;
			typename PolynomialType::PolyType h_dependentPart;

			//Choose Terms from Polynome denpending on dependency on choosen Variable
			for (polynomialIt = inPut.begin(); polynomialIt != inPut.end(); polynomialIt++)
			{
				if (polynomialIt->has(*selectedVariable))
				{
					//divide dependent terms by choosen Variable
					
					polynomialIt->divide(*selectedVariable, tempTerm);
					h_dependentPart.addTerm( tempTerm );
				}
				else 
				{
					h_independentPart.addTerm( *polynomialIt );
				}
			}

			//If Dependent Polynome contains Variables - continue with recursive Horner
			if ( !h_dependentPart.isNumber() )
			{
				#ifdef DEBUG_HORNER
					std::cout << __func__ << "    DEP->new Horner " << h_dependentPart << std::endl;
				#endif
				std::shared_ptr<MultivariateHorner<PolynomialType, Strat>> new_dependent (new MultivariateHorner< PolynomialType, Strat >(std::move(h_dependentPart),s,map));
				setDependent(new_dependent);
				mConst_dependent = constant_zero<CoeffType>::get();	
			}
			
			//Dependent Polynome is a Constant (Number) - save to memberVar
			else 
			{
				removeDependent();
				mConst_dependent = h_dependentPart.constantPart();
			}			

			//If independent Polynome contains Variables - continue with recursive Horner
			if ( !h_independentPart.isNumber() )
			{
				#ifdef DEBUG_HORNER
					std::cout << __func__ << "    INDEP->new Horner " << h_independentPart << std::endl;
				#endif
				std::shared_ptr<MultivariateHorner<PolynomialType, Strat>> new_independent ( new MultivariateHorner< PolynomialType, Strat >(std::move(h_independentPart),s,map));
				setIndependent(new_independent);
				mConst_independent = constant_zero<CoeffType>::get();
			}
			//Independent Polynome is a Constant (Number) - save to memberVar
			else
			{
				removeIndepenent();
				mConst_independent =  h_independentPart.constantPart();
			}					
		} 
					
		//If there are no Variables in the polynome
		else 
		{		
			#ifdef DEBUG_HORNER
				std::cout << __func__ << " [CONSTANT TERM w/o Variables]  " << inPut << std::endl;
			#endif
			removeIndepenent();
			removeDependent();
			mConst_independent = inPut.constantPart();
			this->setVariable( Variable::NO_VARIABLE );
		}
	}



/**
	 * Streaming operator for multivariate HornerSchemes.
	 * @param os Output stream.
	 * @param rhs HornerScheme.
	 * @return `os`.
	 */
template< typename PolynomialType, Strategy Strat > 
std::ostream& operator<<(std::ostream& os,const MultivariateHorner<PolynomialType, Strat>& mvH)
{
	if (mvH.getDependent() && mvH.getIndependent())
	{
		if (mvH.getExponent() != 1)
		{
			return (os << mvH.getVariable() <<"^"<< mvH.getExponent() << " * (" << *mvH.getDependent() << ") + " << *mvH.getIndependent());		
		}
		else
		{
			return (os << mvH.getVariable() << " * (" << *mvH.getDependent() << ") + " << *mvH.getIndependent());			
		}
		
	}
	else if (mvH.getDependent() && !mvH.getIndependent())
	{
		if (mvH.getIndepConstant() == 0)
		{
			if (mvH.getExponent() != 1)
			{
				return (os << mvH.getVariable() <<"^"<< mvH.getExponent() << " * (" << *mvH.getDependent() << ")" );			
			}
			else
			{
				return (os << mvH.getVariable() << " * (" << *mvH.getDependent() << ")" );				
			}
		}
		else
		{
			if (mvH.getExponent() != 1)
			{
				return (os << mvH.getVariable() <<"^"<< mvH.getExponent() << " * (" << *mvH.getDependent() << ") + " << mvH.getIndepConstant());	
			}
			else
			{
				return (os << mvH.getVariable() << " * (" << *mvH.getDependent() << ") + " << mvH.getIndepConstant());		
			}
		}	
	}
	else if (!mvH.getDependent() && mvH.getIndependent())
	{
		if (mvH.getDepConstant() == 1)
		{	
			if (mvH.getExponent() != 1)
			{
				return (os << mvH.getVariable() <<"^"<< mvH.getExponent() << " + " << *mvH.getIndependent() );			
			}
			else
			{
				return (os << mvH.getVariable() << " + " << *mvH.getIndependent() );		
			}			
		}
		else
		{	
			if (mvH.getExponent() != 1)
			{
				return (os << mvH.getDepConstant() << mvH.getVariable() <<"^"<< mvH.getExponent() << " + " << *mvH.getIndependent() );		
			}
			else
			{
				return (os << mvH.getDepConstant() << mvH.getVariable() <<" + " << *mvH.getIndependent() );			
			}
		}
	}	
	else //if (mvH.getDependent() == NULL && mvH.getIndependent() == NULL)
	{
		if (mvH.getVariable() == Variable::NO_VARIABLE)
		{
			return (os << mvH.getIndepConstant());
		}

		else if (mvH.getIndepConstant() == 0)
		{
			if (mvH.getDepConstant() == 1)
			{
				if (mvH.getExponent() != 1)
				{
					return (os << mvH.getVariable() <<"^"<< mvH.getExponent() );	
				}
				else
				{
					return (os << mvH.getVariable());	
				}
				
			}
			else
			{
				if (mvH.getExponent() != 1)
				{
					return (os << mvH.getDepConstant() << mvH.getVariable() <<"^"<< mvH.getExponent());	
				}
				else
				{
					return (os << mvH.getDepConstant() << mvH.getVariable());	
				}
			}
		}
		else 
		{
			if (mvH.getExponent() != 1)
			{
				return (os << mvH.getDepConstant() << mvH.getVariable() <<"^"<< mvH.getExponent() << " + " << mvH.getIndepConstant());	
			}
			else
			{
				return (os << mvH.getDepConstant() << mvH.getVariable() <<" + " << mvH.getIndepConstant());	
			}

		}
	}
}

template<typename PolynomialType, Strategy Strat>
std::shared_ptr<MultivariateHorner<PolynomialType, Strat>> simplify( std::shared_ptr<MultivariateHorner<PolynomialType, Strat>> mvH)
{		

	#ifdef DEBUG_HORNER
		std::cout << __func__ << " Polynome: " << *mvH << std::endl;
	#endif

	if (mvH->getDependent() && (mvH->getDependent()->getDependent() || mvH->getDependent()->getDepConstant() != 0) && !mvH->getDependent()->getIndependent()  && mvH->getDependent()->getIndepConstant() == 0 )
	{
		
		if (mvH->getVariable() == mvH->getDependent()->getVariable())
		{
			mvH->setExponent (mvH->getExponent() + mvH->getDependent()->getExponent()) ;

			if (mvH->getDependent()->getDependent())
			{	
				mvH->setDependent( simplify( mvH->getDependent()->getDependent()) );	
			} 
			else if (mvH->getDependent()->getDepConstant() != 0)
			{
				mvH->setDepConstant(mvH->getDependent()->getDepConstant() );
				mvH->removeDependent();
			}

			if (mvH->getIndependent())
			{	
				mvH->setIndependent( simplify( mvH->getIndependent() ));
			}
	
			return ( simplify(mvH));	
		}
	}

	else if (!mvH->getDependent() && mvH->getIndependent())
	{		
		mvH->setIndependent(simplify ( mvH->getIndependent() ));
		mvH->removeDependent();
		return mvH;
	}

	else if (mvH->getDependent()  && !mvH->getIndependent() )
	{
		mvH->setDependent( simplify ( mvH->getDependent()));
		mvH->removeIndepenent();
		return mvH;
	}
	
	else if (mvH->getDependent() && mvH->getIndependent() )
	{	
		mvH->setDependent( simplify( mvH->getDependent()));
		mvH->setIndependent( simplify ( mvH->getIndependent()));
		return mvH;
	}
	
	return(mvH);
}

}//namespace carl