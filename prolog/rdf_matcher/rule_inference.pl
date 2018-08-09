
:- module(rule_inference,
          [equivalent/2,
           ground_rules/1]).

:- use_module(library(rdf_matcher)).


erule(eq_from_match/5).
erule(eq_from_shared_xref/3).


postulated_rule(equivalent(A,B),Body,P-Rest2) :-
                erule(P/Ar),
                Ar2 is Ar+2,
                functor(Body,P,Ar2),
                Body =.. [P,A,B|Rest],
                subset_nd(Rest,Rest2),
                debug(rdf_matcher,'PRule: ~w ~q',[P,Rest2]).



subset_nd([],[]).
subset_nd([H|T],[H|T2]) :-
        subset_nd(T,T2).
subset_nd([_H|T],['__ANY__'|T2]) :-
        subset_nd(T,T2).




ground_rules(Rules) :-
        debug(rdf_matcher,'Grounding rules',[]),
        setof(FreeVars,H^B^(postulated_rule(H,B,FreeVars),B),AllFreeVars),
        length(AllFreeVars, Num),
        debug(rdf_matcher,'Finding ground rules;  freevars = ~d',[Num]),
        setof((H:-B),BFree^(postulated_rule(H,B,BFree),member(BFree,AllFreeVars)),Rules).

