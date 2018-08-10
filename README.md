# rdf matcher

## Command Line Usage

(see below for use via Docker)

To run these examples, first:

`cd tests/data`


All matches, as TSV:

`rdfmatch  -f tsv -i basic.ttl match`

Contrain to classes with `http://example.org/x/` prefix:

`rdfmatch -p x  -f tsv -i basic.ttl match`

New matches:

`rdfmatch -p x  -f tsv -i basic.ttl new_match`

Cache index in a temp folder:

`rdfmatch -X tmp -p x  -f tsv -i basic.ttl new_match`

Learn match rule probabilities:

`rdfmatch  -p x   -X tmp -T -f tsv -v -i tests/data/basic.ttl -i tests/data/equivs.ttl  learn`

Note this requires a training set - `owl:equivalentClass` is used for this

## Cacheing

The library index_util is used for caching. Optionally this can go to disk

## Docker

Use `./run.sh` instead of `rdfmatch`. E.g.

`./run.sh  -p x  -f tsv -i tests/data/basic.ttl match`

