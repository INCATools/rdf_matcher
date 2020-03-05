#!/bin/sh
docker run -p 9055:9055 -e PORT=9055 -v $PWD:/work -w /work  --rm -ti cmungall/rdf_matcher swipl -G0  -p library=/tools/prolog rdfmatch "$@"

