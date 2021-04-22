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
show_tp(N1,N2,MPred,[Min,Max]) :-
        Opts=[],
        format('T: ~w <-> ~w~n',[N1,N2]),
        assert(fail_test(N1,N2,Pred)),
        forall(term_pair_subsumer_match(N1,N2,Pred,Score,Opts),
               (   format('  M: ~w ~w~n',[Pred,Score]),
                   (   MPred=Pred,
                       Score >= Min,
                       Score =< Max
                   ->  format('     **PASS**~n'),
                       retractall(fail_test(N1,N2,Pred))
                   ;   true))),
        assertion( \+ fail_test(N1,N2,Pred) ).

show_cp(X1,X2,MPred,[Min,Max]) :-
        Opts=[min_score(-8)],
        rdf_global_id(X1,C1),
        rdf_global_id(X2,C2),
        format('TC: ~w <-> ~w~n',[C1,C2]),
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
        show_tp('forebrain development', 'brain development',subClassOf, [4,5]),
        show_tp('forebrain interneuron', 'brain neuron',subClassOf, [3,4]),
        show_tp('brain interneuron', 'forebrain neuron',subClassOf, [0,1]),
        show_tp('brain interneuron', 'forebrain neuron',superClassOf, [0,1]),
        show_tp('H(2)O + thiamine triphosphate = H(+) + phosphate + thiamine diphosphate', 'H2O + thiamine triphosphate = H(+) + phosphate + thiamine diphosphate', equivalentTo,[6,15]),
        true.

test(ont) :-
        show_tp('a hand', 'a manus',equivalentTo, [2,2]),
        show_tp('manus development', 'hand development',equivalentTo, [4,5]),
        show_tp('MANUS DEVELOPMENT', 'hand development',equivalentTo, [4,5]),
        show_tp('development manus', 'hand development',equivalentTo, [4,5]),
        %show_tp('development of manus', 'hand development',equivalentTo, [4,5]),
        show_tp('bone development', 'tissue development',subClassOf, [4,5]),
        true.
% x:bone_of_head  BONE, HEAD      equivalentTo    z:hindlimb      hindlimb        Lexical x       z       rdf_matcher     0.8017474370922646      .


test(cp) :-
        show_cp(x:bone_element,y:bone,equivalentTo,[3,4]),
        show_cp(x:heart,y:heart,equivalentTo,[2,4]),
        show_cp(a:bone_of_foot,x:foot,equivalentTo,[-8,0]),
        G1=basic_annot( 'http://example.org/a/bone_of_foot',_,_,_),
        forall(G1,writeln(G1)),               
        G=tr_annot( 'http://example.org/a/bone_of_foot',_,_,_,_,_),
        forall(G,writeln(G)),               
        assertion( \+ tr_annot( 'http://example.org/a/bone_of_foot',_,foot,_,_,_)),
        %show_cp(x:bone_of_head,z:hindlimb,equivalentTo,[0,0]),
        true.

        



:- end_tests(subsumer_match).
    
