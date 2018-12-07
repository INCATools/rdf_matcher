/**

  given rules or facts with heads such as "eq_from_shared_xref(A,B,X,P1,P2)"

  find all partially grounded rules such as
  
  equivalent(A,B) :- eq_from_shared_xref(A,B,'UBERON',z,x)

  A and B are always unground, but remaining arguments may be a mix
  
  
*/

:- module(rule_inference,
          [equivalent/2,
           ground_rules/1,
           ground_rules/2]).

:- use_module(library(rdf_matcher)).


% assume first two args are class pair
erule(eq_from_match/7).
erule(eq_from_shared_xref/5).


postulated_rule(equivalent(A,B),Body,P-Rest2,Opts) :-
                erule(P/Ar),
                functor(Body,P,Ar),
                Body =.. [P,A,B|Rest],
                (   option(ontA(OntA),Opts)
                ->  reverse(Rest,[_,OntA|_])
                ;   true),
                (   option(ontB(OntB),Opts)
                ->  reverse(Rest,[OntB|_])
                ;   true),
                debug(rdf_matcher,'Candidate Rule: ~q :- ~q',[eq(A,B),Body]),        
                subset_nd(Rest,Rest2),
                debug(rdf_matcher,'PRule: ~w Body=~q // tmpl= ~q',[P,Body,Rest2]).



subset_nd([],[]).
subset_nd([H|T],[H|T2]) :-
        subset_nd(T,T2).
subset_nd([H|T],['__ANY__'|T2]) :-
        var(H),
        subset_nd(T,T2).



ground_rules(Rules) :-
        ground_rules(Rules,[]).
ground_rules(Rules, Opts) :-
        debug(rdf_matcher,'Grounding rules',[]),
        setof(FreeVars,H^B^(postulated_rule(H,B,FreeVars,Opts),B),AllFreeVars),
        length(AllFreeVars, Num),
        debug(rdf_matcher,'Finding ground rules;  freevars = ~d',[Num]),
        setof((H:-B),BFree^(postulated_rule(H,B,BFree,Opts),member(BFree,AllFreeVars)),Rules).

/*
genrule(equivalent(A,B),Body,Opts) :-
                erule(P/Ar),
                functor(Body,P,Ar),
                Body =.. [P,A,B|Rest],
                (   option(ontA(OntA),Opts)
                ->  reverse(Rest,[_,OntA|_])
                ;   true),
                (   option(ontB(OntB),Opts)
                ->  reverse(Rest,[OntB|_])
                ;   true),
                debug(rdf_matcher,'Candidate Rule:~q',[Body]).

ground_rule(Rule) :-
        postulated_rule(H,B,FreeVars,Opts),
        B,
*/        
