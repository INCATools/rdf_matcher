:- use_module(library(rdf_matcher)).
:- use_module(library(rdf_matcher/rule_inference)).
:- use_module(library(semweb/rdf11)).
:- use_module(library(semweb/rdf_turtle)).
:- use_module(library(semweb/rdf_turtle_write)).

:- rdf_register_prefix(a,'http://example.org/a/').
:- rdf_register_prefix(x,'http://example.org/x/').
:- rdf_register_prefix(y,'http://example.org/y/').
:- rdf_register_prefix(z,'http://example.org/z/').

:- debug(subsumer).

:- begin_tests(subsumer_match,
               [setup(load_test_file_and_index),
                cleanup(rdf_retractall(_,_,_,_))]).

load_test_file_and_index :-
        rdf_load('tests/data/basic.ttl'),
        index_pairs.

showall(G) :-
        forall(G,
               format('~q.~n',[G])).


:- dynamic fail_test/3.
:- dynamic max_score/1.
show_tp(N1,N2,MPred,[Min,Max]) :-
        show_tp(N1,N2,MPred,[Min,Max],_).
show_tp(N1,N2,MPred,[Min,Max],MaxScore) :-
        Opts=[min_weight(-8)],
        format('T: ~w <-> ~w // ~w~n',[N1,N2,MPred]),
        assert(fail_test(N1,N2,Pred)),
        nb_setval(max_score,-99),
        forall(term_pair_subsumer_match(N1,N2,Pred,Score,Opts),
               (   format('  M: ~w ~w~n',[Pred,Score]),
                   (   MPred=Pred,
                       Score >= Min,
                       Score =< Max
                   ->  format('     **PASS**~n'),
                       (   nb_getval(max_score,CurrMax),
                           (   Score > CurrMax
                           ->  nb_setval(max_score,Score)
                           ;   true)),
                       retractall(fail_test(N1,N2,Pred))
                       ;   true))),
        nb_getval(max_score,MaxScore),
        assertion( \+ fail_test(N1,N2,Pred) ).

show_cp(X1,X2,MPred,[Min,Max]) :-
        Opts=[min_weight(-8)],
        rdf_global_id(X1,C1),
        rdf_global_id(X2,C2),
        format('TC: ~w <-> ~w  Expect: ~w~n',[C1,C2,MPred]),
        assert(fail_test(C1,C2,Pred)),
        forall(class_pair_subsumer_match(C1,C2,Pred,Score,Opts),
               (   format('  M: ~w ~w~n',[Pred,Score]),
                   (   MPred=Pred,
                       Score >= Min,
                       Score =< Max
                   ->  format('     **PASS**~n'),
                       retractall(fail_test(C1,C2,Pred))
                   ;   true))),
        assertion( \+ fail_test(C1,C2,Pred) ).


test(simple) :-
        show_tp('abc', 'abc',equivalentTo, [1,2]),
        show_tp('abc', 'ABC',equivalentTo, [1,2]),
        show_tp('abc def', 'abc def',equivalentTo, [2,3]),
        show_tp('abc def', 'ab ef', subClassOf, [1,2]),
        show_tp('ab ef', 'abc def', superClassOf, [1,2]),
        show_tp('forebrain development', 'brain development',subClassOf, [3,5]),
        show_tp('forebrain interneuron', 'brain neuron',subClassOf, [3,4]),  % both terms more specific on left
        show_tp('brain interneuron', 'forebrain neuron',subClassOf, [0,1]),   % mixture of more specific and more general
        assertion( \+ term_pair_subsumer_match('brain interneuron', 'forebrain neuron',superClassOf,_,[])),
        show_tp('H(2)O + thiamine triphosphate = H(+) + phosphate + thiamine diphosphate', 'H2O + thiamine triphosphate = H(+) + phosphate + thiamine diphosphate', equivalentTo,[6,15]),
        true.

test(ont) :-
        show_tp('a hand', 'a manus',equivalentTo, [2,2]),
        show_tp('manus development', 'hand development',equivalentTo, [4,5]),
        show_tp('MANUS DEVELOPMENT', 'hand development',equivalentTo, [4,5]),
        show_tp('development manus', 'hand development',equivalentTo, [4,5]),
        show_tp('development of manus', 'hand development',equivalentTo, [3,4]),
        show_tp('bone development', 'tissue development',subClassOf, [4,5]),
        show_tp('bone development', 'tissue development',equivalentTo, [1,2]),
        true.


test(cp) :-
        show_cp(x:bone_element,y:bone,equivalentTo,[3,4]),
        show_cp(x:heart,y:heart,equivalentTo,[2,4]),
        show_cp(a:bone_of_foot,x:foot,equivalentTo,[-8,0]),
        assertion( \+ tr_annot( 'http://example.org/a/bone_of_foot',_,foot,_,_,_)),
        assertion(\+ class_pair_subsumer_match(x:bone_of_head,z:hindlimb,equivalentTo,_,[])),
        true.

test(comparative) :-
        show_tp('a', 'a',equivalentTo, [0,0.5]),
        show_tp('alpha beta gamma', 'alpha beta gamma',equivalentTo, [4,8], S3m), % 3 matches
        show_tp('alpha beta', 'beta gamma',equivalentTo, [-4,0], S3x1),
        assertion(S3m > S3x1),        
        show_tp('alpha beta gamma delta epsilon omega', 'alpha beta gamma delta epsilon omega',equivalentTo, [8,16], S6m),
        assertion(S6m > S3m),        
        show_tp('alpha beta gamma delta epsilon omega', 'alpha beta gamma delta epsilon omexx',equivalentTo, [8,16], S6x5),
        assertion(S6m > S6x5),
        show_tp('xxalpha beta gamma delta epsilon omega', 'alpha beta gamma delta epsilon omexx',equivalentTo, [6,10], S6x4),
        assertion(S6x5 > S6x4),
        show_tp('alpha beta gamma delta', 'alpha beta gamma',subClassOf, [4,8], S3mSub),
        show_tp('alpha beta gamma delta epsilon', 'alpha beta gamma',subClassOf, [4,8], S3mSub2),
        %assertion(S3mSub > S3mSub2),
        show_tp('alpha beta gamma delta epsilon', 'alpha beta gamma phi',subClassOf, [2,4], S3mSub3),
        assertion(S3mSub > S3mSub3),
        assertion(S3mSub2 > S3mSub3),
        true.
        



:- end_tests(subsumer_match).
    
