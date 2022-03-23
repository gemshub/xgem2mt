#!/bin/bash
# Installing dependencies needed to build xGEM2MT for linux 18.04

if [ "$(uname)" == "Darwin" ]; then

    # Do under Mac OS X platform
    EXTN=dylib

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then

    EXTN=so
fi

sudo rm -f /usr/local/lib/libGEMS3K.$EXTN
sudo rm -f /usr/local/lib/libGEM2MT.$EXTN

threads=1
BRANCH_GEMS=xgm2mt-use
BuildType=Release

# GEMS3K
# if no GEMS3K installed in /usr/local/lib/libGEMS3K.so (/usr/local/include/GEMS3K)
# if no GEM2MT installed in /usr/local/lib/libGEM2MT.so (/usr/local/include/nodearray-gem)
test -f /usr/local/lib/libGEM2MT.$EXTN || {

        # Building GEMS3K library
        mkdir -p ~/code && \
                cd ~/code && \
                git clone https://bitbucket.org/gems4/gems3k.git  -b $BRANCH_GEMS && \
                cd gems3k && \
                mkdir -p build && \
                cd build && \
                cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DBUILD_GEM2MT=ON -DCMAKE_BUILD_TYPE=$BuildType && \
                make -j $threads && \
                sudo make install

        # Removing generated build files
        #cd ~ && \
        #         rm -rf ~/code
}

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
   sudo ldconfig
fi

