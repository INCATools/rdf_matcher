FROM swipl:7.5.12

RUN apt-get update && apt-get install make

# Run the image as a non-root user
RUN useradd -m -s /bin/sh myuser
USER myuser
WORKDIR /home/myuser

ADD . $HOME

EXPOSE ${PORT}

## RUN swipl -g "getenv('HOME',Home),atom_concat('file://',Home,Path),Opts=[interactive(false)],pack_install(Path,Opts),halt"
RUN swipl -g "Opts=[interactive(false)],pack_install(index_util,Opts),pack_install(sparqlprog,Opts),halt"
##RUN swipl -g "Opts=[interactive(false)],pack_install(index_util,Opts),halt"

CMD swipl -p library=prolog ./bin/rdfmatch -h
