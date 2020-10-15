[![DOI](https://zenodo.org/badge/13996/cmungall/rdf_matcher.svg)](https://zenodo.org/badge/latestdoi/13996/cmungall/rdf_matcher)

# rdf matcher

## Command Line Usage

(see below for use via Docker)

To run these examples, first:

`cd tests/data`


All matches, as TSV:

`rdfmatch  -f tsv -i basic.ttl match`

Also save triples

`rdfmatch  -f tsv -i basic.ttl -G matches.ttl match`

As above, use skos rather than owl as vocab

`rdfmatch  -f tsv -i basic.ttl -G matches.ttl --predicate skos:exactMatch match`

Contrain to classes with `http://example.org/x/` prefix:

`rdfmatch -p x  -f tsv -i basic.ttl match`

New matches:

`rdfmatch -p x  -f tsv -i basic.ttl new_match`

Cache index in a temp folder:

`rdfmatch -X tmp -p x  -f tsv -i basic.ttl new_match`

Learn match rule probabilities:

`rdfmatch  -p x   -X tmp -T -f tsv -v -i tests/data/basic.ttl -i tests/data/equivs.ttl  learn > probrules.pro`

Note this requires a training set - `owl:equivalentClass` is used for this

Once you you created the probabilistic rules, these can be applied:

`rdfmatch -f tsv -i tests/data/basic.ttl classify probrules.pro`

Generate mappings in SSSOM format:

`rdfmatch  -w prolog/rdf_matcher/obo_weights.pl -T -i tests/data/basic.ttl match`




## Cacheing

The library index_util is used for caching. Optionally this can go to disk

## Docker

Use `./run.sh` instead of `rdfmatch`. E.g.

`./run.sh  -p x  -f tsv -i tests/data/basic.ttl match`

Note that files must be in current directory or below, [run.sh](run.sh) automaps the current dir to the container.
