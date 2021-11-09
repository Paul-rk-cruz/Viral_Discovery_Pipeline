FROM ubuntu:latest
MAINTAINER Paul RK Cruz
LABEL version="1.0"
LABEL description="Viral Discovery Pipeline v1.0 - Docker Dependencies "

# INSTALL SYSTEM DEPENDENCIES
RUN apt-get clean
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install build-essential -y
RUN apt install git -y
RUN apt-get install -y perl
RUN apt-get install -y libexpat1-dev
RUN apt-get install -qy git
RUN apt-get install -qy locales
RUN apt-get install -qy nano
RUN apt-get install -qy tmux
RUN apt-get install -qy wget
RUN apt-get install -qy python3
RUN apt-get install -qy python3-psycopg2
RUN apt-get install -qy python3-pystache
RUN apt-get install -qy python3-yaml
RUN apt-get install -qy curl
RUN apt-get install -qy cpanminus

# INSTALL HOMEBREW & PIPELINE DEPENDENCIES
RUN git clone https://github.com/Homebrew/brew ~/.linuxbrew/Homebrew \
&& mkdir ~/.linuxbrew/bin \
&& ln -s ../Homebrew/bin/brew ~/.linuxbrew/bin \
&& eval $(~/.linuxbrew/bin/brew shellenv) \
&& brew install diamond \
&& brew install spades \
&& brew --version 
 
# CLEANUP DOCKER FILE SYSTEM
RUN apt-get -qy autoremove

CMD ["/bin/bash"]