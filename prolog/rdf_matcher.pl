
:- module(rdf_matcher,
          [
           infer_labels/0,
           obj/1,
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

           match_is_inexact/1,

           create_bitmap_index/0,
           atom_bm/3,
           atom_semsim_match/4,
           pair_semsim_match/5,

           class_prop_bm/4,
           bm_tokens/2,
           bm_resnik/4,

           inter_pair_subsumer_match/5,
           class_pair_subsumer_match/5,           
           term_pair_subsumer_match/5,           
           pair_match/4,
           pair_cmatch/4,
           intra_pair_match/4,
           inter_pair_match/4,
           inter_pair_cmatch/4,
           exact_inter_pair_match/11,
           new_inter_pair_match/4,
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


           categorize_match/3,
           
           eq_from_match/7,
           eq_from_shared_xref/5,

           unmatched/1,
           unmatched_in/2,

           remove_non_english_literals/0,
           remove_language_literals_except/1,
           remove_inexact_synonyms/0,
           transfer_annotation_provenance/0,
           inferred_stopword/1,
           inferred_stopword/2,
           inferred_stopword/3,
           inferred_stopword/4,

           synthesize_class/5,
           merge_into/3,
           add_match_to_graph/3,

           entity_category/2
           ]).

:- use_module(library(porter_stem)).
:- use_module(library(isub)).
:- use_module(library(index_util)).
%:- use_module(library(tabling)).
:- use_module(library(semweb/rdf11)).
:- use_module(library(semweb/rdfs)).
:- use_module(library(sparqlprog/owl_util)).
:- use_module(library(sparqlprog/search_util)).

:- use_module(library(settings)).
:- setting(ontology, atom,'','').

:- rdf_register_prefix(skos, 'http://www.w3.org/2004/02/skos/core#').
:- rdf_register_prefix(skosxl, 'http://www.w3.org/2008/05/skos-xl#').
:- rdf_register_prefix(inca, 'https://w3id.org/inca/').
:- rdf_register_prefix(schema, 'http://schema.org/').
:- rdf_register_prefix(biopax3, 'http://www.biopax.org/release/biopax-level3.owl#').


:- multifile user:weight/2.
:- dynamic user:weight/2.

set_ontology(Ont) :-
        debug(rdf_matcher, 'Setting ontology to ~w',[Ont]),
        set_setting(ontology, Ont).


:- rdf_register_prefix(oio,'http://www.geneontology.org/formats/oboInOwl#').
%:- rdf_register_prefix(skos,'http://www.w3.org/2004/02/skos/core#"').

%:- module_transparent index_pairs/0, index_pairs/1.
index_pairs :-
        index_pairs(none).
index_pairs(no_index) :- !.
index_pairs(Path) :-
        %materialize_index(obj(+)),
        Goals=[
               %equivalent(+,+),
               basic_annot(+,+,-,-),
               tr_annot(+,+,+,-,-,-)
               %pair_match(+,+,-,-),
               %inter_pair_match(+,+,-,-)
              ],
        (   Path=none
        ->  maplist(materialize_index,Goals)
        ;   materialize_indexes_to_path(Goals, Path)).



% TODO: make easier to extend/plugin

:- rdf_meta pmap(-,r).

pmap(label, skos:prefLabel).
pmap(label, rdfs:label).
pmap(label, biopax3:standardName).
pmap(related, skos:altLabel).
pmap(related, oio:hasRelatedSynonym).
pmap(related, biopax3:name).
pmap(exact, oio:hasExactSynonym).
pmap(broad, oio:hasBroadSynonym).
pmap(narrow, oio:hasNarrowSynonym).

pmap(related, skos:closeMatch).
pmap(exact, skos:exactMatch).
pmap(broad, skos:broadMatch).
pmap(narrow, skos:narrowMatch).

pmap(distinct, skos:relatedMatch).

pmap(xref, oio:hasDbXref).

%objpmap(exact, skos:exactMatch).

inexact(broad).
inexact(narrow).
inexact(related).

match_is_inexact(info(P-_,_,_)) :- inexact(P),!.
match_is_inexact(info(_-P,_,_)) :- inexact(P),!.
match_is_inexact(info(_,_,stem)) :- !.
        

nonmut(xref).
nonmut(id).
nonmut(uri).


rdf_object_as_atom(X,X) :- rdf_is_iri(X),!.
rdf_object_as_atom(X,A) :- literal_atom(X,A),!.



literal_string(S^^_,S).
literal_string(S@_,S).

literal_atom(L,A) :- literal_string(L,S),!,atom_string(A,S).

opt_literal_atom(L,A) :- literal_atom(L,A), !.
opt_literal_atom(A,A) :- atomic(A).

:- table obj/1.
obj(Obj) :-
        setof(Obj, X^Y^rdf(Obj, X , Y), Objs),
        member(Obj,Objs),
        rdf_is_iri(Obj),
        has_prefix(Obj,Prefix),
        \+ is_builtin_prefix(Prefix),
        \+ skos_xl_instance(Obj).


is_builtin_prefix(rdfs).
is_builtin_prefix(http).

skos_xl_instance(Obj) :-
        rdf(Obj,rdf:type,C),
        rdf_global_id(skosxl:_,C).

%obj(Obj) :-
%        setof(Obj, rdf(Obj,rdf:type,owl:'Class'), Objs),
%        member(Obj,Objs),
%        rdf_is_iri(Obj).


%% basic_annot(?Object, ?AnnotProp, ?Val ,?RdfTripleTerm) is nondet
basic_annot(Obj,P,V) :-
        basic_annot(Obj,P,V,_).
basic_annot(Obj,P,V,T) :-
        pmap(P,P1),
        obj(Obj),
        T=rdf(Obj,P1,Val),
        T,
        rdf_object_as_atom(Val,V).
basic_annot(Obj,id,V,rdf(Obj,dc:identifier,V)) :-
        obj(Obj),
        rdf_global_id(Pre:Post,Obj),
        concat_atom([Pre,Post],:,V).
basic_annot(Obj,uri,Obj,rdf(Obj,schema:url,Obj)) :-
        obj(Obj).
%basic_annot(Obj,P,Val,uri(Val)) :-
%        objpmap(P,P1),
%        obj(Obj),
%        T=rdf(Obj,P1,Val),
%        T.

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


mutate(tokenset,_,V,V2) :-
        tokenset_atom(V,V2).
mutate(stem,_,V,V2) :-
        custom_porter_stem(V,V2).
mutate(downcase,_,V,V2) :-
        downcase_atom(V,V2).

%TODO: this fails in a random fashion
%:- table term_replacement/2.
term_replacement(Pattern,Replacement) :-
        rdf(S,inca:pattern,PatternLit),
        rdf(S,inca:replacement,ReplacementLit),
        literal_atom(PatternLit,Pattern),
        literal_atom(ReplacementLit,Replacement).
xxxterm_replacement(Pattern,'') :-
        inferred_stopword(Pattern),
        atom_chars(Pattern,Cs),
        forall(member(C,Cs),
               is_alpha(C)).
term_replacement(eous,eus). % TODO config
term_replacement('  *',' ').
term_replacement('^ +','').
term_replacement(' +$','').

replace_termsyns(T,S) :-
        nb_setval(term,T),
        forall(term_replacement(Pat,Rep),
               (   nb_getval(term,T1),
                   re_replace(Pat,Rep,T1,T2),
                   nb_setval(term,T2))),
        nb_getval(term,S),
        !.
replace_termsyns(T,_) :-
        throw(error(no_replace(T))).


% TODO: tokenize

custom_porter_stem(T,S) :-
        is_ascii(T),
        replace_termsyns(T,T2),
        porter_stem(T2,S).

tokenset_atom(T,S) :-
        is_ascii(T),
        downcase_atom(T,T1),
        replace_termsyns(T1,T2),
        tokenize_atom(T2,As),
        sort(As,As2),
        include(is_true_token,As2,As3),
        concat_atom(As3,S).

is_true_token(A):-
        sub_atom(A,0,1,_,Letter),
        is_alphanumeric(Letter).

is_alphanumeric(X) :-  X @>= 'a',  X @=< 'z',!.
is_alphanumeric(X) :-  X @>= 'A',  X @=< 'Z',!.
is_alphanumeric(X) :-  X @>= '0',  X @=< '9'.


        

% porter_stem only accepts ISO-Latin 1. To avoid conversion
% (not sure how) we simply do not attempt to stem anything that isn't on
% the ascii subset
is_ascii(T) :-
        atom_codes(T,Codes),
        \+ ((member(C,Codes), C > 255)).

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


%! obj_token(?Obj,?Tok) is det
%
%   maps between an OWL object IRI and a token obtained from tokenizing a label or synonym
%
obj_token(Obj,Tok) :-
        obj_token(Obj,Tok,_,_).
obj_token(Obj,Tok,P,V) :-
        rdf(Obj,rdf:type,_),
        tr_annot(Obj,P,V,_,_,_),
        \+ nonmut(P),
        every_n(objtoken,10000,debug(index,'OT: ~w ~w ~w',[Obj,P,V])),
        concat_atom(Toks,' ',V),
        member(Tok,Toks),
        atom_length(Tok,TokLen),
        % TODO: make this less arbitrary
        % designed to filter out things that are not truly word tokens,
        % or hard to tokenize chemical names. http://purl.obolibrary.org/obo/CHEBI_36054
        TokLen < 25.

every_n(Name,Num,Goal) :-
        (   nb_current(Name,X)
        ->  X2 is X+1
        ;   X2=1),
        nb_setval(Name,X2),
        Mod is X2 mod Num,
        (   Mod=0
        ->  debug(counter,'~w Iteration ~w',[Name,X2]),
            Goal
        ;   true),
        !.

        
:- multifile index_util:fast_index/1.

index_util:fast_index(obj_token/4).
index_util:fast_index(obj_token/2).
index_util:fast_index(token_index/2).
index_util:fast_index(position_ic/2).

:- dynamic position_ic/2.
create_bitmap_index :-
        %assert(token_index_stored('',0)),
        debug(index,'Getting all obj-token pairs ot/4',[]),
        materialize_index(obj_token(+,+,+,+)),
        debug(index,'Getting all obj-token pairs ot/2',[]),
        materialize_index(obj_token(+,+)),
        debug(index,'getting distinct tokens...',[]),
        %setof(Token, Obj^obj_token(Obj,Token),Tokens),
        findall(Token, obj_token(Obj,Token),TokensNU),
        length(TokensNU,NumUsages),
        debug(index,'uniquifying... Usages: ~w',[NumUsages]),
        sort(TokensNU,Tokens),
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

combine_mutfunc(X,X,X) :- !.
combine_mutfunc(null,downcase,null) :- !.
combine_mutfunc(downcase,null,null) :- !.

myopt(Opt,Opts,Default) :-
        Opt =.. [_,A],
        option(Opt,Opts,Default),
        (   var(A)
        ->  A=Default
        ;   true).

cls_lexform(C,N) :-
        tr_annot(C,A,N,_,_,_),
        debug(subsumer,'trann: ~q ~q ~q',[C,A,N]),
        A\=uri,
        A\=xref,
        A\=id.
inter_pair_subsumer_match(C1,C2,Pred,Score,Opts) :-
        class_pair_subsumer_match(C1,C2,Pred,Score,Opts),
        has_prefix(C1,Pfx1),
        has_prefix(C2,Pfx2),
        Pfx1 \= Pfx2.
class_pair_subsumer_match(C1,C2,Pred,Score,Opts) :-
        cls_lexform(C1,N1),
        cls_lexform(C2,N2),
        debug(subsumer,'CP: ~w <-> ~w ( ~w ~w) // "~w" <-> "~w"',[C1,C2,A1,A2,N1,N2]),
        term_pair_subsumer_match(N1,N2,Pred,Score,Opts).

:- table stopword/1.
stopword('') :- !.
stopword(W) :- rdf(_,inca:stopword,Lit),
        literal_atom(Lit,W).

:- table mytok/2.
mytok(N,Toks3) :-
        concat_atom(Toks,' ',N),
        maplist(downcase_atom,Toks,Toks2),
        findall(T,(member(T,Toks2),
                   \+stopword(T)),
                Toks3).
term_pair_subsumer_match(N1,N2,Pred,Score,Opts) :-
        mytok(N1,Toks1),
        mytok(N2,Toks2),
        sort(Toks1,Set1),
        sort(Toks2,Set2),
        tokens_intersections(Set1,Set2,Ixn,Set1U,Set2U),
        length(Ixn,NumIxn),
        debug(subsumer,'~q <-> ~q :: Ixn=~q',[Set1,Set2,Ixn]),
        myopt(min_ixn(MinIxn),Opts,0),
        NumIxn >= MinIxn,
        findall(match(X,X,equivalentTo,W),(member(X,Ixn),identical_match_weight(X,W)),Matches1),
        debug(subsumer,'  Matches1 :: ~q',[Matches1]),
        tokens_subsumers(Set1U,Set2U,Matches2),
        append(Matches1,Matches2,Matches),
        debug(subsumer,'  Matches :: ~q',[Matches]),
        member(Pred,[subClassOf,superClassOf,equivalentTo]),
        \+ \+ member(match(_,_,Pred,_),Matches),
        findall(W,(member(M,Matches),smatch_weight(M,Pred,W)),Weights),
        debug(subsumer,'  Weights :: ~q',[Weights]),
        sum_list(Weights,Score),
        myopt(min_score(MinScore),Opts,0),
        Score > MinScore.


smatch_weight(match(_,_,P,W),P,W) :- !.
smatch_weight(match(_,_,_,W),equivalentTo,W2) :- W2 is W/2, !.
smatch_weight(match(_,_,equivalentTo,W),P,W) :- !.
smatch_weight(match(_,_,_,_),_,-1).
smatch_weight(nomatch(_),_,-2).

identical_match_weight(M,W) :-
        identical_match_weight1(M,W),
        !.
identical_match_weight1(a,0).
identical_match_weight1(to,0).
identical_match_weight1(the,0).
identical_match_weight1(X,W) :-
        atom_length(X,Len),
        W is log(Len + 1).


tokens_intersections(Set1,Set2,Ixn,Set1U,Set2U) :-
        ord_intersection(Set1,Set2,Ixn),
        ord_subtract(Set1,Ixn,Set1U),
        ord_subtract(Set2,Ixn,Set2U).
tokens_subsumers([],Set2,Matches) :-
        !,
        findall(nomatch(X),member(X,Set2),Matches).
tokens_subsumers(Set1,Set2,[Match|Matches]) :-
        select(E1,Set1,Set1R),
        best_subsumer(E1,Set2,Set2R,Match),
        !,
        tokens_subsumers(Set1R,Set2R,Matches).
tokens_subsumers(Set1,Set2,[nomatch(E1)|Matches]) :-
        select(E1,Set1,Set1R),
        tokens_subsumers(Set1R,Set2,Matches).

best_subsumer(E1,Set2,Set2R,Match) :-
        findall(Score-Match,(member(E2,Set2),token_subsumer(E1,E2,Score,Match)),Matches),
        sort_best(Matches,_-Match),
        match(_,E2,_,_) = Match,
        select(E2,Set2,Set2R).

:- table token_subsumer/4.
token_subsumer(E1,E2,S,match(E1,E2,Rel,S)) :-
        token_subsumer_nondir(E1,E2,S,Rel).
token_subsumer(E1,E2,S,match(E1,E2,IRel,S)) :-
        token_subsumer_nondir(E2,E1,S,Rel),
        invert(Rel,IRel).

token_subsumer_nondir(E1,E2,Score,superClassOf) :-
        sub_atom(E2,_,MLen,_,E1),
        atom_length(E2,Len),
        Extra is log(MLen),  % longer matches better
        Score is ((MLen / Len) - 0.5) + Extra.
token_subsumer_nondir(E1,E2,Score,equivalentTo) :-
        isub(E1,E2,true,Score),
        Score > 0.5,
        debug(subsumer,'xxxxxxxxxxISUB ~q <-> ~q = ~w',[E1,E2,Score]).
token_subsumer_nondir(E1,E2,2,equivalentTo) :-
        cls_lexform(C,E1),
        cls_lexform(C,E2).
token_subsumer_nondir(E1,E2,Score,Rel) :-
        cls_lexform(C1,E1),
        cls_lexform(C2,E2),
        node_subsumer_nondir(C1,C2,Score,Rel).

node_subsumer_nondir(C1,C2,1.8,subClassOf) :-
        rdf(C1,rdfs:subClassOf,C2).
node_subsumer_nondir(C1,C2,0.95,subClassOf) :-
        rdfs_subclass_of(C1,C2),
        \+ rdf(C1,rdfs:subClassOf,C2).





% TODO: ontology subsumption

invert(equivalentTo,equivalentTo).
invert(subClassOf,superClassOf).
invert(superClassOf,subClassOf).



sort_best(L,Best) :-
        sort(L,LS),
        reverse(LS,[Best|_]).


%% pair_match(?Class1, ?Class2, ?SharedVal, Info) is nondet
%
% Info = ?AP1, ?AP2, ?Triple1, ?Triple2, ?MutFunc
:- rdf_meta pair_match(r,r,-,-).
:- rdf_meta pair_match(r,r).
pair_match(C1,C2,V,Info) :-
        Info=info(P1-P2,T1-T2,MutFunc),
        tr_annot(C1,P1,V,T1,MutFunc1,_),
        tr_annot(C2,P2,V,T2,MutFunc2,_),
        combine_mutfunc(MutFunc1,MutFunc2,MutFunc),
        \+ excluded(C1,C2),
        C1\=C2.
pair_match(C1,C2) :- pair_match(C1,C2,_,_).


pair_cmatch(C1,C2,V,Info) :-
        rdf_global_id(C1,C1x),
        rdf_global_id(C2,C2x),
        pair_match(C1x,C2x,V,Info).

:- rdf_meta intra_pair_match(r,r,-,-).
intra_pair_match(C1,C2,V,Info) :-
        pair_match(C1,C2,V,Info),
        C1\=C2.
%has_prefix(C1,Pfx),
%has_prefix(C2,Pfx).
intra_pair_match(null,null,null,null).


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

new_inter_pair_match(C1,C2,V,Info) :-
        inter_pair_match(C1,C2,V,Info),
        \+ equivalent(C1,_),
        \+ equivalent(_,C1),
        \+ equivalent(C2,_),
        \+ equivalent(_,C2).


%:- rdf_meta new_pair_match(r,r,-,-).
%
% an inter_pair_match/4 in which the prefixes
% for which there is no existing match for C1 in prefix(C2)
% and there is no existing match for C2 in prefix(C1)
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


        
exact_inter_pair_match(C,X,CParents,XParents,Conf,Vs,Info,AltCs,AltXs,IgnoredCs,IgnoredXs) :-
        setof(V-Info,
              exact_inter_pair_match1(C,X,CParents,XParents,Conf,V,Info,AltCs,AltXs,IgnoredCs,IgnoredXs),
              Ms),
        setof(V,Info^member(V-Info,Ms),Vs),
        setof(Info,V^member(V-Info,Ms),Infos),
        sformat(Info,'~w',[Infos]).

exact_inter_pair_match1(C,X,CParents,XParents,Conf,V,info(Pred1,Pred2,Func),AltCs,AltXs,IgnoredCs,IgnoredXs) :-
        exact_inter_pair_match(C,X,V,Info),
        Info=info(Pred1-Pred2,_,Func),
        findall(X2,alt_exact_inter_pair_match(C,X,X2),AltXs),
        findall(C2,alt_exact_inter_pair_match(X,C,C2),AltCs),
        findall(X2,alt_inexact_inter_pair_match(C,X,X2),IgnoredXs),
        findall(C2,alt_inexact_inter_pair_match(X,C,C2),IgnoredCs),
        (   AltXs=[],
            AltCs=[]
        ->  (   IgnoredXs=[],
                IgnoredCs=[]
            ->  Conf=high
            ;   Conf=medium)
        ;   Conf=low),
        findall(Parent,entity_parent(C,Parent),CParents),
        findall(Parent,entity_parent(X,Parent),XParents).

exact_inter_pair_match(C,X,V,Info) :-
        inter_pair_match(C,X,V,Info),
        \+ match_is_inexact(Info).

entity_parent(X,Parent) :-
        rdf(X,rdfs:subClassOf,Parent),
        rdf_is_iri(Parent).
entity_parent(X,Parent) :-
        subclass_of_some(X,R,Parent),
        parent_relation(R),
        rdf_is_iri(Parent).
entity_parent(X,Parent) :-
        rdf(X,R,Parent),
        parent_relation(R),
        rdf_is_iri(Parent).
entity_parent(X,Parent) :-
        rdf(X,rdf:type,Parent),
        \+ rdf_global_id(rdf:_,Parent),
        \+ rdf_global_id(rdfs:_,Parent),
        \+ rdf_global_id(owl:_,Parent).


alt_exact_inter_pair_match(C,X,X2) :-
        exact_inter_pair_match(C,X2,_,_),
        X2\=X.
alt_inexact_inter_pair_match(C,X,X2) :-
        inter_pair_match(C,X2,_,Info),
        X2\=X,
        match_is_inexact(Info).



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
        has_prefix(X,P),
        debug(prefixes,'Prefix ~w -> ~w',[X,P]).

guess_uribase(X,U) :-
        rdf(X,_,_),
        defrag(X,U,_).

defrag(X,U,Frag) :-
        concat_atom(Parts,'/',X),
        reverse(Parts,[Frag|Rev]),
        reverse(Rev,Parts2),
        concat_atom(Parts2,'/',U).


has_alternative_obj_match(Sub,Obj,true) :-
        has_prefix(Obj,ObjSrc),
        pair_match(Sub,Obj2,_,_),
        Obj2\=Obj,
        has_prefix(Obj2,ObjSrc),
        !.
has_alternative_obj_match(_Sub,_Obj,false).

has_alternative_sub_match(Sub,Obj,true) :-
        has_prefix(Sub,SubSrc),
        pair_match(Sub2,Obj,_,_),
        Sub2\=Sub,
        has_prefix(Sub2,SubSrc),
        !.
has_alternative_sub_match(_Sub,_Obj,false).


categorize_match(Sub,Obj,Cat) :-
        has_alternative_sub_match(Sub,Obj,AltSub),
        has_alternative_obj_match(Sub,Obj,AltObj),
        member(AltSub+AltObj = Cat,
               [false+false = unique,
                false+true = one_to_many,
                true+false = many_to_one,
                true+true = many_to_many]).



declare_additional_prefixes :-
        rdf(X,'http://www.w3.org/ns/shacl#prefix',^^(Prefix1,_)),
        rdf(X,'http://www.w3.org/ns/shacl#namespace', ^^(NS1,_)),
        atom_string(Prefix,Prefix1),
        atom_string(NS,NS1),
        rdf_register_prefix(Prefix, NS),
        debug(rdf_matcher,'Registered ~w ~w',[Prefix,NS]),
        fail.
declare_additional_prefixes.


inject_prefixes :-
        declare_additional_prefixes,
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

:- table obj_has_prefix/2.
obj_has_prefix(C,P) :- obj(C),has_prefix(C,P).

:- table has_prefix/2.
has_prefix(C,P) :- atomic(C), rdf_global_id(P:_, C).
has_prefix(C,P) :- atomic(C), \+ rdf_global_id(_:_, C), concat_atom([P,_],:,C), P\=http.
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

:- table object_token/2.
:- table object_token/3.
object_token(Obj,T) :-
        object_token(Obj,T,_).
object_token(Obj,T,Prefix) :-
        obj(Obj),
        has_prefix(Obj,Prefix),
        rdf(Obj,P1,Val),
        pmap(_,P1),
        literal_atom(Val,A),
        tokenize_atom(A,Ts),
        member(T,Ts).

:- table token/1.
token(T) :-
        setof(T,O^object_token(O,T),Ts),
        member(T,Ts).
:- table prefix_token/2.
prefix_token(P,T) :-
        setof(T,O^object_token(O,T,P),Ts),
        member(T,Ts).

prefix_token_count(Prefix,T,Num) :-
        aggregate(count,O,object_token(O,T,Prefix),Num).

:- table inferred_stopword/1.
inferred_stopword(T) :-
        inferred_stopword(T,_).
inferred_stopword(T,Prefix) :-
        inferred_stopword(T,Prefix,_).
inferred_stopword(T,Prefix,Prop) :-
        inferred_stopword(T,Prefix,Prop,_).
inferred_stopword(T,Prefix,Prop,Min) :-
        (   var(Min)
        ->  Min = 0.75
        ;   Min =< 1.0),
        aggregate(count,O,obj(O),NumObjsTotal),
        debug(rdf_matcher,'Total objs: ~w',[NumObjsTotal]),
        % arbitrary... vocabs like agrovoc are too big for tabling
        NumObjsTotal < 70000,
        debug(rdf_matcher,'Counting objs for: ~w',[Prefix]),
        aggregate(count,O,obj_has_prefix(O,Prefix),NumObjs),
        debug(rdf_matcher,'Prefix: ~w Objs: ~w',[Prefix,NumObjs]),
        prefix_token_count(Prefix,T,NumWithT),
        debug(stopword,'  T: ~w N: ~w',[T,NumWithT]),
        Prop is NumWithT/NumObjs,
        (   Prop > Min
        ;   common_word(T)).


common_word(the).
             
% todo - figure out why this is so slow
transfer_annotation_provenance :-
        findall(S-XL,triple_axiom_annotation(S,_P,_O,oio:hasDbXref,XL),Pairs),
        forall((member(S-XL,Pairs),
                literal_atom(XL,X),
                concat_atom([Prefix,Local],':',X),
                used_prefix(Prefix),
                rdf_global_id(Prefix:Local,X_URI)),
               (   rdf_assert(S,skos:exactMatch,X_URI),
                   debug(rdf_matcher,'Made skos triple ~w ->  ~w // from annotation',[S,X_URI]))).


remove_inexact_synonyms :-
        T=rdf(_,P,_),
        findall(T,
                (   inexact(PN),
                    pmap(PN,P),
                    T),
                Ts),
        forall(member(rdf(S,P,O),Ts),
               (   debug(rdf_matcher,'Removing: ~w ~w ~w',[S,P,O]),
                   rdf_retractall(S,P,O))).

remove_non_english_literals :-
        remove_language_literals_except(en).
remove_language_literals_except(PreserveLang) :-
        findall(rdf(S,P,Lit),
                (   rdf(S,P,Lit),
                    Lit= _ @ Lang,
                    Lang\=PreserveLang),
                Ts),
        forall(member(rdf(S,P,O),Ts),
               (   debug(rdf_matcher_debug,'Removing: ~w ~w ~w',[S,P,O]),
                   rdf_retractall(S,P,O))).

mpred(skos:exactMatch).
mpred(skos:closeMatch).
mpred(owl:equivalentClass).
mpred(owl:equivalentProperty).
mpred(owl:disjointWith).
asserted_match(X,Y) :- mpred(P),rdf(X,P,Y).
asserted_match(X,Y) :- mpred(P),rdf(Y,P,X).

synthesize_class(Label,Prefix,C,TargetGraph,Opts) :-
        rdf_assert(TargetGraph,rdf:type,owl:'Ontology',TargetGraph),
        rdf_assert(oio:source,rdf:type,owl:'AnnotationProperty',TargetGraph),
        rdf_global_id(Prefix:Label,C),
        rdf_assert(C,rdf:type,owl:'Class'),
        rdf_assert(C,rdfs:label,Label@en),
        setof(X,lsearch(Label,X),Sources),
        merge_into(C, Sources, TargetGraph, Opts).
        

%! merge_all_equivalents(+TargetPrefix, +TargetGraph, ?Opts) is det.
merge_all_equivalents(TargetPrefix, TargetGraph, Opts) :-
        setof(Tgt,obj_has_prefix(Tgt,TargetPrefix),Tgts),
        forall(member(Tgt,Tgts),
               merge_into(Tgt, TargetGraph, Opts)).
        
%! merge_into(+Target:iri, +TargetGraph:iri, +Opts:list) is det.
%! merge_into(+Target:iri, +Sources:list, +TargetGraph:iri, +Opts:list) is det.
merge_into(Target, TargetGraph, Opts) :-
        (   member(match(true),Opts)
        ->  find_and_store_matches(Target,0,matches,Opts)
        ;   true),
        findall(X,asserted_match(Target,X),Sources),
        merge_into(Target, Sources, TargetGraph, Opts).
merge_into(Target, Sources, TargetGraph, Opts) :-
        % bring over axioms from Source, rewiring if necessary
        forall(member(S,Sources),
               forall(rdf(S,P,O),
                      rdf_assert_mapped(Target, P, O, TargetGraph, S, Opts))),
        % bring over axioms about target, unaltered
        forall(rdf(Target,P,O),
               rdf_assert(Target, P, O, TargetGraph)),
        % inject skos exact matches
        forall(member(S,Sources),
               rdf_assert(Target,skos:exactMatch,S, TargetGraph)).

rdf_assert_mapped(S, P, O, TG, SourceEntity, Opts) :-
        map_predicate(P, Px, Opts),
        map_iri(O, Ox, P, Opts),
        (   S=Ox
        ->  true
        ;   rdf_assert_annotated(S, Px, Ox, TG, oio:source, SourceEntity)).

:- rdf_meta map_predicate(r,r,o).

map_predicate(rdfs:label,oio:hasExactSynonym,_) :- !.  % preserve cardinality
map_predicate(P,Px,Opts) :-
        owl_equivalent_property_asserted_symm(P,Px),
        has_prefix(Px,Prefix),
        member(prefix(Prefix),Opts),
        !.
map_predicate(P,P,_).

map_iri(A,A,P,_) :- mpred(Px),rdf_global_id(Px,P),!.  % do not map if object of a mapping
map_iri(A,A,_,_) :- \+ rdf_is_iri(A).
map_iri(A,Ax,_,Opts) :-
        owl_equivalent_class_asserted_symm(A,Ax),
        has_prefix(Ax,Prefix),
        debug(rdf_matcher,'Eq: ~w == ~w ;; checking prefix ~w',[A,Ax,Prefix]),
        member(prefix(Prefix),Opts),
        debug(rdf_matcher,'Mapped ~w -> ~w',[A,Ax]),
        !.
map_iri(A,A,_,_).


% det
find_and_store_matches(X,Dist,Graph,Opts) :-
        option(search_distance(MaxDist),Opts,1),
        Dist < MaxDist,
        match_goal_template(X,Y,G,Opts),
        Dist1 is Dist+1,
        forall(G,(add_match_to_graph(G,Graph,[match_predicate(skos:exactMatch)|Opts]),
                  find_and_store_matches(Y,Dist1,Graph,Opts))).
find_and_store_matches(_,_,_,_).


match_goal_template(X,Y,G,Opts) :-
        member(quickmatch(true),Opts),
        !,
        G=quick_match(X,Y,_).
match_goal_template(X,Y,G,_Opts) :- G = pair_match(X,Y,_,_).

quick_match(X,Y,Info) :-
        rdf(X,XP,XV),
        pmap(XPN,XP),
        literal_atom(XV,XVA),
        {icase(YV,XVA)},
        rdf(Y,YP,YV),
        pmap(YPN,YP),
        Info=m(XPN,YPN).

% for a given successful Goal, serialize to results rdf graph
% NOTE: due to a bug in earlier owlapi a new version should be used see https://github.com/owlcs/owlapi/issues/875
add_match_to_graph(G,Graph,Opts) :-
        debug(rdf_matcher,'Adding match: ~q to ~q // Opts=~w',[G,Graph,Opts]),
        G =.. [_,Sub,Obj|Args],
        option(match_predicate(Pred),Opts,owl:equivalentClass),
        %(   Pred=owl:equivalentClass,
        %    rdf(Sub,rdf:type,Type),
        %    rdf(Obj,rdf:type,Type),
        %    fix_equiv_pred(Type,Pred1)
        %->  true
        %;   Pred1=Pred),
        rdf_global_id(Pred,Pred2),
        rdf_assert(Sub,Pred2,Obj,Graph),
        rdf_create_bnode(Axiom),
        rdf_assert(Axiom,rdf:type,owl:'Axiom',Graph),
        rdf_assert(Axiom,owl:annotatedSource,Sub,Graph),
        rdf_assert(Axiom,owl:annotatedTarget,Obj,Graph),
        rdf_assert(Axiom,owl:annotatedProperty,Pred2,Graph),
        assert_label(Axiom,Graph,inca:sourceLabel,Sub),
        assert_label(Axiom,Graph,inca:targetLabel,Obj),
        (   member(info(_,rdf(_,SourcePred,_)-rdf(_,TargetPred,_),ProcessStep),Args),
            has_prefix(Sub,SubPrefix),
            has_prefix(Obj,ObjPrefix),
            rdf_global_id(SubPrefix:'',SubPrefixURL),
            rdf_global_id(ObjPrefix:'',ObjPrefixURL)
        ->  rdf_assert(Axiom,inca:sourcePredicate,SourcePred,Graph),
            rdf_assert(Axiom,inca:targetPredicate,TargetPred,Graph),
            rdf_assert(Axiom,inca:sourcePrefix,SubPrefixURL,Graph),
            rdf_assert(Axiom,inca:targetPrefix,ObjPrefixURL,Graph),
            rdf_assert(Axiom,inca:processingStep,ProcessStep,Graph)
        ;   true),
        sformat(Info,'Info: ~w',[Args]),
        rdf_canonical_literal(Info,InfoLit),
        rdf_assert(Axiom,rdfs:comment,InfoLit,Graph),
        % this is redundantly asserted each time
        rdf_assert(Graph,rdf:type,owl:'Ontology',Graph).


assert_label(Axiom,Graph,AP,X) :-
        rdf(X,rdfs:label,Label),
        !,
        rdf_assert(AP,rdfs:type,owl:'AnnotationProperty',Graph),
        rdf_assert(Axiom,AP,Label,Graph).
assert_label(_,_,_,_).

rdf_assert_annotated(Sub,Pred,Obj,Graph,P,V) :-
        rdf_create_bnode(Axiom),
        rdf_assert(Axiom,rdf:type,owl:'Axiom',Graph),
        rdf_assert(Axiom,owl:annotatedSource,Sub,Graph),
        rdf_assert(Axiom,owl:annotatedTarget,Obj,Graph),
        rdf_assert(Axiom,owl:annotatedProperty,Pred,Graph),
        rdf_assert(Sub,P,V,Graph).

infer_labels :-
        rdf_subject(C),
        \+ rdf(C,rdfs:label,_),
        infer_label(C,N),
        rdf_assert(C,rdfs:label,N@en),
        fail.
infer_labels.

infer_label(URI,N) :-
        rdf_subject(URI),
        concat_atom([_Base,Frag],'#',URI),
        expand_camel_case(Frag,N1),
        concat_atom(L,'_',N1),
        concat_atom(L,' ',N),
        !.
expand_camel_case(In,Out) :-
        atom_chars(In,InC),
        expand_camel_case_chars(InC,OutC,true),
        atom_chars(Out,OutC).
expand_camel_case_chars([],[],_).
expand_camel_case_chars([H|T],[' ',H|T2],false) :-
        is_camel_hump(H),
        !,
        expand_camel_case_chars(T,T2,false).
expand_camel_case_chars([H|T],[H|T2],_) :-
        !,
        expand_camel_case_chars(T,T2,false).

is_camel_hump(H) :-
        H @>= 'A',
        H @=< 'Z'.
is_camel_hump(H) :-
        H @>= '0',
        H @=< '9'.

foo__('0').

:- table entity_category/2.

%! entity_category(?E,?C) is nondet.
%! entity_category(+E,?C) is semidet.
%
%  if there are multiple categories, |-concatenated
%
%  note this doesn't usually happen as we try and find the most relevant category
%
entity_category(E,C) :-
        setof(C1,(obj(E),entity_category1(E,C1)),Cs),
        !,
        maplist(category_label,Cs,C2s),
        concat_atom(C2s,'|',C).
entity_category(_,'').

category_label(C,N) :-
        atom(C),
        rdf(C,rdfs:label,Label),
        literal_atom(Label,N),
        !.
category_label(C,N) :-
        literal_atom(C,N),
        !.
category_label(C,N) :-
        atom_concat('https://w3id.org/biolink/vocab/', N, C),
        !.
category_label(C,C).

has_skos_top_concept(E,C) :-
        rdf_reachable(E,skos:broader,C1),
        rdf(C1,skos:topConceptOf,_),
        rdf(C1,skos:prefLabel,C).


entity_category1(E,C) :-
        setof(C1,has_skos_top_concept(E,C1),Cs),
        !,
        member(C,Cs).
entity_category1(E,C) :-
        rdf(E,skos:inScheme,C),
        !.
entity_category1(E,C) :-
        % assumes OWL layering of RDF
        rdf(E,rdf:type,owl:'NamedIndividual'),
        rdf(E,rdf:type,C),
        \+ rdf_global_id(owl:_,C),
        !.
%entity_category1(E,C) :-
%        % assumes OWL layering of RDF
%        rdf(E,rdf:type,owl:'NamedIndividual'),
%        rdf(E,rdf:type,T),
%        entity_category1(T,C),
%        !.
entity_category1(E,C) :-
        rdfs_subclass_of(E,P),
        parent_category(P,C,Priority),
        \+ ((rdfs_subclass_of(E,P1),
             rdfs_subclass_of(P1,P),
             P1\=P,
             parent_category(P,C2,Priority2),
             Priority2 < Priority,
             C2\=C)),
        !.
entity_category1(E,C) :-
        rdf(E,oio:hasOBONamespace,C),
        !.
entity_category1(E,C) :-
        % C is a root class
        rdfs_subclass_of(E,C),
        rdf_is_iri(C),
        \+ rdf_global_id(owl:_,C),
        \+ ((rdf(C,rdfs:subClassOf,P),
             \+ rdf_global_id('http://schema.org/Thing',P),
             \+ rdf_global_id(owl:'Thing',P))),
        !.
entity_category1(_,thing) :- !.

parent_category(E,C,1) :-
        category_property(P),
        rdf(E,P,C).
parent_category(E,C,2) :-
        rdf(E,rdf:type,C),
        \+ rdf_global_id(owl:_,C),
        \+ rdf_global_id(skos:_,C),
        \+ rdf_global_id(rdfs:_,C).

category_property('https://w3id.org/biolink/vocab/category').
category_property('http://dbpedia.org/ontology/category').

% make configurable by having equiv axiom
parent_relation(X) :- rdf(X,skos:exactMatch,rdfs:subClassOf).
parent_relation('http://purl.obolibrary.org/obo/gaz#located_in').
parent_relation('http://purl.obolibrary.org/obo/RO_0001025').
parent_relation('http://purl.obolibrary.org/obo/BFO_0000050').


