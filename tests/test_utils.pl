:- module(test_utils,
          []).

:- use_module(library(semweb/rdf11)).

register_prefixes :-
        rdf_register_prefix(x,'http://example.org/x/'),
        rdf_register_prefix(y,'http://example.org/y/'),
        rdf_register_prefix(z,'http://example.org/z/').

:- initialization(register_prefixes, now).

    
