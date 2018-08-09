:- module(rdf_matcher,
          [
           index_pairs/0,
           index_pairs/1,
           obj_has_prefix/2,
           equivalent/2,
           set_ontology/1,
           
           tr_annot/6,
           has_prefix/2,
           used_prefix/1,
           inject_prefixes/0,

           
           pair_match/4,
           pair_cmatch/4,
           inter_pair_match/4,
           inter_pair_cmatch/4,
           new_pair_match/4,
           new_pair_cmatch/4,
           new_unique_pair_match/4,
           new_unique_pair_cmatch/4,
           new_ambiguous_pair_match/6,
           new_ambiguous_pair_cmatch/6,

           transitive_match/2,
           transitive_match_set/1,
           transitive_match_set_member/2,

           eq_from_match/7,
           eq_from_shared_xref/5
           ]).

:- use_module(library(porter_stem)).
:- use_module(library(index_util)).
:- use_module(library(tabling)).
:- use_module(library(semweb/rdf11)).

:- use_module(library(settings)).
:- setting(ontology, atom,'','').

set_ontology(Ont) :-
        debug(rdf_matcher, 'Setting ontology to ~w',[Ont]),
        set_setting(ontology, Ont).


:- rdf_register_prefix(oio,'http://www.geneontology.org/formats/oboInOwl#').
%:- rdf_register_prefix(skos,'http://www.w3.org/2004/02/skos/core#"').

%:- module_transparent index_pairs/0, index_pairs/1.
index_pairs :-
        index_pairs(none).
index_pairs(Path) :-
        materialize_index(obj(+)),
        Goals=[
               %equivalent(+,+),
               basic_annot(+,+,-,-),
               tr_annot(+,+,+,-,-,-),
               pair_match(+,+,-,-),
               inter_pair_match(+,+,-,-)
              ],
        (   Path=none
        ->  maplist(materialize_index,Goals)
        ;   materialize_indexes_to_path(Goals, Path)).



% TODO: make easier to extend/plugin

:- rdf_meta pmap(-,r).

pmap(label, rdfs:label).
pmap(related, oio:hasRelatedSynonym).
pmap(exact, oio:hasExactSynonym).
pmap(broad, oio:hasBroadSynonym).
pmap(narrow, oio:hasNarrowSynonym).

pmap(related, skos:closeMatch).
pmap(exact, skos:exactMatch).
pmap(broad, skos:broadMatch).
pmap(narrow, skos:narrowMatch).

pmap(xref, oio:hasDbXref).

nonmut(xref).
nonmut(id).


literal_string(S^^_,S).
literal_string(S@_,S).

literal_atom(L,A) :- literal_string(L,S),atom_string(A,S).

opt_literal_atom(L,A) :- literal_atom(L,A), !.
opt_literal_atom(A,A) :- atomic(A).

obj(Obj) :-
        setof(Obj, rdf(Obj,rdf:type,owl:'Class'), Objs),
        member(Obj,Objs),
        rdf_is_iri(Obj).


%% basic_annot(?Object, ?AnnotProp, ?Val ,?RdfTripleTerm) is nondet
basic_annot(Obj,P,V) :-
        basic_annot(Obj,P,V,_).
basic_annot(Obj,P,V,T) :-
        pmap(P,P1),
        obj(Obj),
        T=rdf(Obj,P1,Lit),
        T,
        literal_atom(Lit,V).
basic_annot(Obj,id,V,id(V)) :-
        obj(Obj),
        rdf_global_id(Pre:Post,Obj),
        concat_atom([Pre,Post],:,V).
basic_annot(Obj,uri,Obj,uri(Obj)) :-
        obj(Obj).

%% tr_annot(?Object, ?AnnotProp, ?MutVal ,?RdfTripleTerm, ?MutFunc, ?OrigVal) is nondet
tr_annot(Obj,P,V2,T,F,V) :-
        basic_annot(Obj,P,V,T),
        \+ nonmut(P),
        mutate(F,P,V,V2).
tr_annot(Obj,P,V,T,null,V) :-
        nonmut(P),
        basic_annot(Obj,P,V,T).


mutate(stem,_,V,V2) :-
        custom_porter_stem(V,V2).
mutate(downcase,_,V,V2) :-
        downcase_atom(V,V2).

% TODO: tokenize
custom_porter_stem(T,S) :-
        atom_concat(Obj,eous,T),
        atom_concat(Obj,eus,T2),
        !,
        porter_stem(T2,S).
custom_porter_stem(T,S) :-
        porter_stem(T,S).

excluded(C1,_) :-
        setting(ontology,P),
        P\='',
        \+ has_prefix(C1,P).


%% pair_match(?Class1, ?Class2, ?SharedVal, Info) is nondet
%
% Info = ?AP1, ?AP2, ?Triple1, ?Triple2, ?MutFunc
:- rdf_meta pair_match(r,r,-,-).
:- rdf_meta pair_match(r,r).
pair_match(C1,C2,V,Info) :-
        Info=info(P1-P2,T1-T2,MutFunc),
        tr_annot(C1,P1,V,T1,MutFunc,_),
        tr_annot(C2,P2,V,T2,MutFunc,_),
        \+ excluded(C1,C2),
        C1\=C2.
pair_match(C1,C2) :- pair_match(C1,C2,_,_).

pair_cmatch(C1,C2,V,Info) :-
        rdf_global_id(C1,C1x),
        rdf_global_id(C2,C2x),
        pair_match(C1x,C2x,V,Info).


:- rdf_meta inter_pair_match(r,r,-,-).
inter_pair_match(C1,C2,V,Info) :-
        pair_match(C1,C2,V,Info),
        has_prefix(C1,Pfx1),
        has_prefix(C2,Pfx2),
        Pfx1 \= Pfx2.

inter_pair_cmatch(C1,C2,V,Info) :-
        rdf_global_id(C1,C1x),
        rdf_global_id(C2,C2x),
        inter_pair_match(C1x,C2x,V,Info).


:- rdf_meta new_pair_match(r,r,-,-).
new_pair_match(C1,C2,V,Info) :-
        inter_pair_match(C1,C2,V,Info),
        has_prefix(C1,Pfx1),
        has_prefix(C2,Pfx2),
        unmatched_in(C1, Pfx2),
        unmatched_in(C2, Pfx1).

new_pair_cmatch(C1,C2,V,Info) :-
        rdf_global_id(C1,C1x),
        rdf_global_id(C2,C2x),
        new_pair_match(C1x,C2x,V,Info).


:- rdf_meta new_unique_pair_match(r,r,-,-).
new_unique_pair_match(C1,C2,V,Info) :-
        new_pair_match(C1,C2,V,Info),
        has_prefix(C1,Pfx1),
        has_prefix(C2,Pfx2),
        \+ ((pair_match(C1,AltC2,_,_)),
            AltC2\=C2,
            has_prefix(AltC2,Pfx2)),
        \+ ((pair_match(AltC1,C2,_,_)),
            AltC1\=C1,
            has_prefix(AltC1,Pfx1)).

new_unique_pair_cmatch(C1,C2,V,Info) :-
        rdf_global_id(C1,C1x),
        rdf_global_id(C2,C2x),
        new_unique_pair_match(C1x,C2x,V,Info).


new_ambiguous_pair_match(C1,C2,AltC1,AltC2,V,Info) :-
        new_pair_match(C1,C2,V,Info),
        has_prefix(C1,Pfx1),
        has_prefix(C2,Pfx2),
        (   pair_match(C1,AltC2,_,_),
            AltC2\=C2,
            has_prefix(AltC2,Pfx2),
            AltC1='-'
        ;   pair_match(AltC1,C2,_,_),
            AltC1\=C1,
            has_prefix(AltC1,Pfx1),
            AltC2='-').

new_ambiguous_pair_cmatch(C1,C2,AltC1,AltC2,V,Info) :-
        rdf_global_id(C1,C1x),
        rdf_global_id(C2,C2x),
        new_ambiguous_pair_match(C1x,C2x,AltC1,AltC2,V,Info).

:- table transitive_match/2.
transitive_match(A,B) :-
        new_unique_pair_match(A,Z,_,_),
        transitive_match(Z,B).
transitive_match(A,B) :-
        new_unique_pair_match(A,B,_,_).

transitive_match_set(Bs) :-
        obj(A),
        setof(B,transitive_match(A,B),Bs),
        Bs=[_,_,_|_],
        \+ ((member(Z,Bs),
            equivalent(Z,_))).

transitive_match_set_member(X,M) :-
        transitive_match_set(Set),
        member(M,Set),
        Set=[X|_].





%! eq_from_match(?ClsA, ?ClsB, ?SynPredA, ?SynPredB, ?MutateOp, ?OntA, ?OntB) is nondet
%
% 
eq_from_match(A,B,APred,BPred,Mut,OntA,OntB) :-
        inter_pair_match(A,B,_,info(APred-BPred, _, Mut)),
        has_prefix(A,OntA),
        has_prefix(B,OntB).
eq_from_shared_xref(A,B,OntX,OntA,OntB) :-
        inter_pair_match(A,B,X,info(xref-xref, _, null)),
        has_prefix(A,OntA),
        has_prefix(B,OntB),
        has_prefix(X,OntX).

used_prefixes(Ps) :-
        setof(P,used_prefix1(P),Ps).
used_prefix(P) :-
        setof(P,used_prefix1(P),Ps),
        member(P,Ps).
used_prefix1(P) :-
        rdf(X,rdf:type,_),
        has_prefix(X,P).

guess_uribase(X,U) :-
        rdf(X,_,_),
        defrag(X,U,_).

defrag(X,U,Frag) :-
        concat_atom(Parts,'/',X),
        reverse(Parts,[Frag|Rev]),
        reverse(Rev,Parts2),
        concat_atom(Parts2,'/',U).
        

inject_prefixes :-
        (   used_prefixes(Ps),
            debug(rdf_matcher,'Found these prefiexes: ~w',[Ps]),
            Ps\=[_,_|_]
        ->  force_inject_prefixes
        ;   true).
force_inject_prefixes :-
        debug(rdf_matcher,'Injecting prefixes',[]),
        setof(U,X^guess_uribase(X,U),Us),
        debug(rdf_matcher,'Bases = ~w',[Us]),
        Us=[_,_|_],
        forall((member(U,Us),defrag(U,_,Prefix)),
               (   debug(rdf_matcher,'Register: ~w : ~w',[Prefix,U]),
                   rdf_register_prefix(Prefix, U))),
        !.
force_inject_prefixes :-
        throw(error(cannot_guess_prefixes)).

        
        
        
equivalent(C1,C2) :- rdf(C1,owl:equivalentClass,C2).
equivalent(C1,C2) :- rdf(C2,owl:equivalentClass,C1).

obj_has_prefix(C,P) :- obj(C),has_prefix(C,P).

has_prefix(C,P) :- atomic(C), rdf_global_id(P:_, C).
has_prefix(C,P) :- atomic(C), \+ rdf_global_id(_:_, C), concat_atom([P,_],:,C).
has_prefix(P:_,P).

%% unmatched(?Cls) is nondet.
%
% true if Cls has no equivalent class
unmatched(C) :-
        obj(C),
        \+ equivalent(C,_).

%% unmatched_in(?Cls, ?ExtPrefix) is nondet.
%
% true if Cls has no equivalent class with prefix ExtPrefix
unmatched_in(C, ExtPrefix) :-
        obj(C),
        \+ ((equivalent(C,C2),
             has_prefix(C2, ExtPrefix))).

