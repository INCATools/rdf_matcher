FROM swipl:7.5.12
MAINTAINER Chris Mungall <cjmungall@lbl.gov>

RUN apt-get update && apt-get install make
ADD . /tools
WORKDIR /tools
RUN swipl -g "Opts=[interactive(false)],pack_install(index_util,Opts),pack_install(sparqlprog,Opts),halt"
ENV PATH "/tools/bin:$PATH"

CMD swipl -p library=prolog ./bin/rdfmatch -h
