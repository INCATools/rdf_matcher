:- module(rdf_matcher,
          [
           index_pairs/0,
           index_pairs/1,
           obj_has_prefix/2,
           equivalent/2,
           set_ontology/1,

           filter_mapped_classes/0,

           mutate/4,
           basic_annot/4,
           tr_annot/6,
           has_prefix/2,
           used_prefix/1,
           inject_prefixes/0,

           create_bitmap_index/0,
           atom_bm/3,
           atom_semsim_match/4,
           pair_semsim_match/5,

           class_prop_bm/4,
           bm_tokens/2,
           bm_resnik/4,
           
           pair_match/4,
           pair_cmatch/4,
           inter_pair_match/4,
           inter_pair_cmatch/4,
           new_pair_match/4,
           new_pair_cmatch/4,
           rightnew_pair_match/4,
           rightnew_pair_cmatch/4,
           new_unique_pair_match/4,
           new_unique_pair_cmatch/4,
           new_ambiguous_pair_match/6,
           new_ambiguous_pair_cmatch/6,

           tri_match/4,
           new_unique_match_triad_nc/3,
           
           transitive_unique_match/2,
           transitive_unique_match_set/1,
           transitive_unique_match_set_member/2,

           transitive_new_match/2,
           transitive_new_match_set_pair/3,

           eq_from_match/7,
           eq_from_shared_xref/5,

           unmatched/1,
           unmatched_in/2
           ]).

:- use_module(library(porter_stem)).
:- use_module(library(index_util)).
:- use_module(library(tabling)).
:- use_module(library(semweb/rdf11)).

:- use_module(library(settings)).
:- setting(ontology, atom,'','').

:- rdf_register_prefix(skos, 'http://www.w3.org/2004/02/skos/core#').

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

%! tr_annot(?Object, ?AnnotProp, ?MutVal ,?RdfTripleTerm, ?MutFunc, ?OrigVal) is nondet
%
%    AnnotProp = id | ...
%
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

:- dynamic token_index_stored/2.
token_index(Tok,Ix) :-
        (   token_index_stored(Tok,Ix)
        *->  true
        ;   nonvar(Tok),
            gensym('',IxAtom),
            atom_number(IxAtom,Ix),
            assert(token_index_stored(Tok,Ix))).

obj_token(Obj,Tok) :-
        obj_token(Obj,Tok,_,_).
obj_token(Obj,Tok,P,V) :-
        rdf(Obj,rdf:type,_),
        tr_annot(Obj,P,V,_,_,_),
        every_n(objtoken,10000,debug(index,'OT: ~w ~w ~w',[Obj,P,V])),
        concat_atom(Toks,' ',V),
        member(Tok,Toks),
        atom_length(Tok,TokLen),
        % TODO: make this less arbitrary
        % designed to filter out things that are not truly word tokens,
        % or hard to tokenize chemical names. http://purl.obolibrary.org/obo/CHEBI_36054
        TokLen < 50.

every_n(Name,Num,Goal) :-
        (   nb_current(Name,X)
        ->  X2 is X+1
        ;   X2=1),
        nb_setval(Name,X2),
        Mod is X2 mod Num,
        (   Mod=0
        ->  debug(counter,'~w Iteration ~w',[Name,X2]),
            Goal
        ;   true).

        

        

:- dynamic position_ic/2.
create_bitmap_index :-
        %assert(token_index_stored('',0)),
        debug(index,'Getting all obj-token pairs',[]),
        materialize_index(obj_token(+,+,+,+)),
        materialize_index(obj_token(+,+)),
        debug(index,'getting distinct tokens...',[]),
        setof(Token, Obj^obj_token(Obj,Token),Tokens),
        debug(index,'aggregating...',[]),
        maplist([Token,Count-Token] >> aggregate(count,Obj,obj_token(Obj,Token),Count), Tokens, CTPairs1),
        debug(index,'sorting...',[]),
        sort(CTPairs1,CTPairs),
        %setof(Count-Token, aggregate(count,Obj,obj_token(Obj,Token),Count), CTPairs),
        debug(index,'getting all counts...',[]),
        findall(Count,member(Count-_,CTPairs),Counts),
        sumlist(Counts,Total),
        debug(index,'Total objs: ~w',[Total]),
        maplist([_-Token,Token-Ix]>>token_index(Token,Ix),CTPairs,_TIPairs),
        %debug(bm,'Pairs=~w',[TIPairs]),
        debug(index,'Materializing token_index',[]),
        materialize_index(token_index(+,+)),
        debug(index,'Calculating ICs',[]),
        forall((member(Count-Token,CTPairs),token_index(Token,Pos)),
               (   IC is -log(Count/Total)/log(2),
                   assert(position_ic(Pos,IC)))),
        debug(index,'Materializing IC index',[]),
        materialize_index(position_ic(+,+)),
        debug(index,'Indexing classes',[]),
        forall(rdf(C,rdf:type,owl:'Class'),index_class(C)),
        debug(index,'Bitmap index complete',[]).
        

:- dynamic class_prop_bm/4.
index_class(C) :-
        setof(Tok,obj_token(C,Tok,P,V),Toks),
        tokens_bm(Toks,BM,_),
        assert(class_prop_bm(C,P,V,BM)),
        every_n(ixc,1000,debug(index,'Class: ~w',[C])),
        fail.
index_class(_).

atom_bm(A,BM,U) :-
        concat_atom(Toks,' ',A),
        tokens_bm(Toks,BM,U).

tokens_bm(Toks,BM,LenU) :-
        maplist([In,N]>>(token_index(In,Ix) -> N is 2**Ix ; N=0),Toks,Nums),
        findall(Tok,(member(Tok,Toks),\+token_index(Tok,_)),UToks),
        length(UToks,LenU),
        sumlist(Nums,BM).

bm_simJ(A,B,S) :-
        I is A /\ B,
        CI is popcount(I),
        (   CI=0
        ->  S=0
        ;   U is A \/ B,
            CU is popcount(U),
            S is CI/CU).

% as above, with unmatched
bm_simJ(A,B,Unmatched,S) :-
        I is A /\ B,
        (   I=0
        ->  S=0
        ;   U is A \/ B,
            CI is popcount(I),
            CU is popcount(U) + Unmatched,
            S is CI/CU).

% fuzzy match of how much A is subsumed by B
%  e.g. a b c subsumed_by a c
%  e.g. a b c d partially subsumed_by a e
bm_subsumed_by_simJ(A,B,S) :-
        I is A /\ B,
        (   I=0
        ->  S=0
        ;   CI is popcount(I),
            CU is popcount(B),
            S is CI/CU).

bm_resnik(A,B,Unmatched,S) :-
        I is A /\ B,
        CI is popcount(I),
        (   CI=0
        ->  S=0
        ;   U is A \/ B,
            bm_sum_ic(I,TI),
            bm_sum_ic(U,TU),
            % default IC for unmatched
            Penalty is Unmatched * 5,
            S is TI/(TU+Penalty)).

bm_subsumed_by_resnik(A,B,S) :-
        I is A /\ B,
        (   I=0
        ->  S=0
        ;   bm_sum_ic(I,TI),
            bm_sum_ic(B,TU),
            S is TI/TU).



atom_semsim_match(A,Cls,S,Method) :-
        atom_semsim_match(A,Cls,_,_,S,Method).
atom_semsim_match(A,Cls,P,V,S,simj) :-
        atom_bm(A,ABM,AU),
        class_prop_bm(Cls,P,V,CBM),
        bm_simJ(ABM,CBM,AU,S).
atom_semsim_match(A,Cls,P,V,S,subsumed_by_simj) :-
        atom_bm(A,ABM,AU),
        class_prop_bm(Cls,P,V,CBM),
        bm_subsumed_by_simj(ABM,CBM,AU,S).
atom_semsim_match(A,Cls,P,V,S,icratio) :-
        atom_bm(A,ABM,AU),
        class_prop_bm(Cls,P,V,CBM),
        bm_resnik(ABM,CBM,AU,S).
atom_semsim_match(A,Cls,P,V,S,subsumed_by_icratio) :-
        atom_bm(A,ABM,_),
        class_prop_bm(Cls,P,V,CBM),
        bm_subsumed_by_resnik(ABM,CBM,S).

pair_semsim_match(icratio,C1,C2,Info,S) :-
        Info=info(P1-P2,V1-V2,u),
        class_prop_bm(C1,P1,V1,BM1),
        class_prop_bm(C2,P2,V2,BM2),
        bm_resnik(BM1,BM2,0,S).


bm_sum_ic(BM,SumIC) :-
        bm_positions(BM,Posns),
        maplist([Pos,IC]>>position_ic(Pos,IC),Posns,ICs),
        sumlist(ICs,SumIC).

bm_tokens(BM, Toks) :-
        bm_positions(BM, Posns),
        maplist([Pos,Tok]>>token_index(Tok,Pos), Posns, Toks).


%% bm_positions(+AV:int,?AL:list)
% True if AV is an integer bit vector with the attributes in AL set
bm_positions(AV,AL) :-
        bm_positions(AV,AL,65536).

bm_positions(AV,AL,Window) :-
        Mask is 2**Window -1,
        bm_positions(AV,ALx,0,Window,Mask),
        flatten(ALx,AL).

%% bm_positions(+AV:int,?AL:list,+Pos,+Window,+Mask) is det
% Mask must = Window^2 -1 (not checked)
% shifts AV down Window bits at a time. If there are any bits in the window,
% use bm_positions_lo/2 to get the attribute list from this window.
% note resulting list must be flattened.
% todo: difference list impl?
bm_positions(0,[],_,_,_) :- !.
bm_positions(AV,AL,Pos,Window,Mask) :-
        !,
        NextBit is AV /\ Mask,
        AVShift is AV >> Window,
        NextPos is Pos+Window,
        (   NextBit=0
        ->  bm_positions(AVShift,AL,NextPos,Window,Mask)
        ;   bm_positions_lo(NextBit,ALNew,Pos),
            AL=[ALNew|AL2],
            bm_positions(AVShift,AL2,NextPos,Window,Mask)).

:- table bm_positions_lo/2.

% as bm_positions/2, but checks one bit at a time
bm_positions_lo(AV,AL) :-
        bm_positions_lo(AV,AL,0).

bm_positions_lo(0,[],_) :- !.
bm_positions_lo(AV,AL,Pos) :-
        NextBit is AV /\ 1,
        AVShift is AV >> 1,
        NextPos is Pos+1,
        (   NextBit=1
        ->  AL=[Pos|AL2]
        ;   AL=AL2),
        !,
        bm_positions_lo(AVShift,AL2,NextPos).

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
inter_pair_match(null,null,null,null).

inter_pair_cmatch(C1,C2,V,Info) :-
        rdf_global_id(C1,C1x),
        rdf_global_id(C2,C2x),
        inter_pair_match(C1x,C2x,V,Info).


%:- rdf_meta new_pair_match(r,r,-,-).
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

% satisfied if C1 and C2 match
% AND C2 has no other match in the same namespace as C1
:- rdf_meta rightnew_pair_match(r,r,-,-).
rightnew_pair_match(C1,C2,V,Info) :-
        inter_pair_match(C1,C2,V,Info),
        has_prefix(C1,Pfx1),
        unmatched_in(C2, Pfx1).

rightnew_pair_cmatch(C1,C2,V,Info) :-
        rdf_global_id(C1,C1x),
        rdf_global_id(C2,C2x),
        rightnew_pair_match(C1x,C2x,V,Info).


% satisfied if C1 and C2 match
% AND C2 has no other match in the same namespace as C1
% AND C1 has no other match in the same namespace as C2
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

tri_match(C1,C2,C3,Info) :-
        new_unique_pair_match(C1,C2,_,_),
        new_unique_pair_match(C2,C3,_,_),
        C3\=C1,
        (   equivalent(C1,C3)
        ->  Info=agree
        ;   (   new_unique_pair_match(C1,C3,_,_)
            ->  Info=all_match
            ;   Info=mismatch)).

        

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

/*
  UNIQUE MATCH CLUSTERS
*/

%:- table transitive_unique_match/2.
transitive_unique_match(A,B) :-
        new_unique_pair_match(A,Z,_,_),
        transitive_unique_match(Z,B).
transitive_unique_match(A,B) :-
        new_unique_pair_match(A,B,_,_).


new_unique_match_triad_nc(A,B,C) :-
        new_unique_pair_match(A,B,_,_),
        A @< B,
        \+ equivalent(A,_),
        \+ equivalent(B,_),
        new_unique_pair_match(B,C,_,_),
        B @< C,
        \+ equivalent(C,_).




transitive_unique_match_set(Bs) :-
        obj(A),
        recursive_expand(A,transitive_unique_match,Bs),
        %setof(B,transitive_unique_match(A,B),Bs),
        Bs=[_,_,_|_],
        \+ ((member(Z,Bs),
            equivalent(Z,_))).


%! transitive_unique_match_set_member(?X, ?M) is nondet
%
%  X is the reference member of a clique, and M is a member
%
transitive_unique_match_set_member(X,M) :-
        transitive_unique_match_set(Set),
        member(M,Set),
        Set=[X|_].

/*
  NEW MATCH CLUSTERS
*/

%:- table transitive_new_match/2.
transitive_new_match(A,B) :-
        new_pair_match(A,Z,_,_),
        transitive_new_match(Z,B).
transitive_new_match(A,B) :-
        new_pair_match(A,B,_,_).

transitive_new_match_set(Bs) :-
        obj(A),
        recursive_expand(A,transitive_new_match,Bs),
        Bs=[_,_,_|_].
transitive_new_match_set_member(X,M) :-
        transitive_new_match_set(Set),
        member(M,Set),
        Set=[X|_].
transitive_new_match_set_pair(X,A,B) :-
        transitive_new_match_set(Set),
        Set=[X|_],
        member(A,Set),
        member(B,Set),
        A@<B,
        new_pair_match(A,B,_,_).


recursive_expand(A,Pred,RSet) :-
        set_recursive_expand([A],Pred,[],RSet).
set_recursive_expand([],_,RSet,RSet2) :-
        sort(RSet,RSet2).
set_recursive_expand([A|Seeds],Pred,Visited,RSet) :-
        Goal=..[Pred,A,X],
        (   setof(X,(Goal,\+member(X,Visited)),Xs)
        ->  ord_union(Seeds,Xs,Seeds2),
            set_recursive_expand(Seeds2,Pred,[A|Visited],RSet)
        ;   set_recursive_expand(Seeds,Pred,[A|Visited],RSet)).


%! eq_from_match(?ClsA, ?ClsB, ?SynPredA, ?SynPredB, ?MutateOp, ?OntA, ?OntB) is nondet
%
%   true if pair ClsA and ClsB match using the given predicate pair (e.g. exact,exact)
%   after MutateOp performed, and belong to OntA and OntB (prefixes)
eq_from_match(A,B,APred,BPred,Mut,OntA,OntB) :-
        inter_pair_match(A,B,_,info(APred-BPred, _, Mut)),
        has_prefix(A,OntA),
        has_prefix(B,OntB).

%! eq_from_shared_xref(?ClsA, ?ClsB, ?SharedXrefOnt, ?OntA, ?OntB) is nondet
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
            debug(rdf_matcher,'Found these prefixes: ~w',[Ps]),
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

        
filter_mapped_classes :-
        setof(C,C2^equivalent(C,C2),MappedCs),
        forall(member(C,MappedCs),
               rdf_retractall(C,_,_,_)).
        
equivalent(C1,C2) :- rdf(C1,owl:equivalentClass,C2).
equivalent(C1,C2) :- rdf(C2,owl:equivalentClass,C1).
equivalent(C1,C2) :- rdf(C1,skos:exactMatch,C2).
equivalent(C1,C2) :- rdf(C2,skos:exactMatch,C1).

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

