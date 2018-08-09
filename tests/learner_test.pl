:- use_module(library(rdf_matcher)).
:- use_module(library(rdf_matcher/rule_inference)).
:- use_module(library(rule_eval)).
:- use_module(library(semweb/rdf11)).
:- use_module(library(semweb/rdf_turtle)).
:- use_module(library(semweb/rdf_turtle_write)).

:- rdf_register_prefix(x,'http://example.org/x/').
:- rdf_register_prefix(y,'http://example.org/y/').
:- rdf_register_prefix(z,'http://example.org/z/').

:- debug(index).
:- debug(rule_eval).

:- begin_tests(learner,
               [setup(load_test_file_and_index),
                cleanup(rdf_retractall(_,_,_,_))]).

load_test_file_and_index :-
        rdf_load('tests/data/basic.ttl'),
        rdf_load('tests/data/equivs.ttl'),
        index_pairs.




show_rules(Rs) :-
        eval_rules(Rs,Pairs),
        forall(member(P,Pairs),
               show_rule(P)).

show_rule(R-Scores) :-
        wrule(R),
        format('Scores: ~q.~n',[Scores]),
        nl.

wrule(R) :-
        copy_term(R,R2),
        numbervars(R2,0,_,[]),
        format('RULE: ~q.~n',[R2]).

        


test(learn) :-
        ground_rules(Rules),
        %setof(FreeVars,H^B^(postulated_rule(H,B,FreeVars),B),AllFreeVars),
        %setof((H:-B),BFree^(postulated_rule(H,B,BFree),member(BFree,AllFreeVars)),Rules),
        show_rules(Rules),
        nl.

:- end_tests(learner).
    
