:- use_module(library(rdf_matcher)).
:- use_module(library(rdf_matcher/rule_inference)).
:- use_module(library(semweb/rdf11)).
:- use_module(library(semweb/rdf_turtle)).
:- use_module(library(semweb/rdf_turtle_write)).

:- rdf_register_prefix(x,'http://example.org/x/').
:- rdf_register_prefix(y,'http://example.org/y/').
:- rdf_register_prefix(z,'http://example.org/z/').

:- debug(index).

:- begin_tests(rdf_matcher,
               [setup(load_test_file_and_index),
                cleanup(rdf_retractall(_,_,_,_))]).

load_test_file_and_index :-
        rdf_load('tests/data/basic.ttl'),
        index_pairs.

showall(G) :-
        forall(G,
               format('~q.~n',[G])).

test(load) :-
        rdf_global_id(x:lung, XLUNG),
        assertion(tr_annot(XLUNG, label, lung, _, _, _)),
        assertion(tr_annot(XLUNG, xref, 'UBERON:0002048', _, _, _)),
        assertion(tr_annot(XLUNG, id, 'x:lung', _, _, _)),

        rdf_global_id(y:lung, YLUNG),
        assertion(tr_annot(YLUNG, label, lungs, _, downcase, _)),
        assertion(tr_annot(YLUNG, label, lung, _, stem, _)),
        
        %showall(tr_annot(_,_,_,_,_,_)),
        %showall(inter_pair_cmatch(_,_,_,_)),
        %showall(new_unique_pair_cmatch(_,_,_,_)),
        %showall(new_ambiguous_pair_cmatch(_,_,_,_,_,_)),

        % match using URIs
        assertion( inter_pair_match('http://example.org/z/bone_tissue','http://example.org/y/bone', _, _) ),
        
        % match using CURIEs
        assertion( inter_pair_cmatch(x:organ, y:organ, _, _) ),

        % match via broad synonym, no stemming
        assertion( inter_pair_cmatch(x:bone_tissue, y:bone, bone, info(broad-label,_,stem)) ),

        % no self-matches
        assertion( \+ pair_cmatch(x:eye, x:eye, _, _)),

        % symmetry
        assertion( (pair_cmatch(x:eye, y:eye, _, _), pair_cmatch(y:eye, x:eye, _, _)) ),
                
        % intra-ontology matches: include with pair_cmatch, exclude with inter_pair_cmatch
        assertion( pair_cmatch(x:bone_tissue, y:bone, _, _)), 

        % 
        assertion( inter_pair_cmatch(x:organ, y:organ, _, info(_,_,stem)) ),

        assertion( new_ambiguous_pair_cmatch(x:bone_tissue, y:bone, _, _, _, _)),

        forall( (inter_pair_match(A,B,_,_), A@<B),
                rdf_assert(A,owl:equivalentClass,B, equivGraph)),

        rdf_save_turtle('tests/data/foo.ttl',[graph(equivGraph),format(ttl)]),
              
        
        nl.

:- end_tests(rdf_matcher).
    
