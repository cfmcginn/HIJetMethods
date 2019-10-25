#!/bin/bash

#### 1. Set top directory of project ####
#Opt for this more complicated solution over PWD so invocation anywhere will still give top dir of HIJetMethods; effectively tho, just=$PWD
internalHIJETMETHODS="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#Let's add some very basic protection in case user has moved the setHIJetMethodsEnv.sh file
testHYDROSTR="HIJetMethods"
tempHIJETMETHODS=${internalHIJETMETHODS%/}
while [[ $tempHIJETMETHODS == *"/"* ]]
do
    tempHIJETMETHODS=${tempHIJETMETHODS#*/}
done

if [[ $testHYDROSTR != $tempHIJETMETHODS ]]
then
    echo "HIJETMETHODS '$internalHIJETMETHODS' DOES NOT TERMINATE AT '$testHYDROSTR'. Please make sure setHIJetMethodsEnv.sh is at top level of HIJetMethods repo. exit 1, HIJETMETHODS unset"
    exit 1    
fi

#We passed check, so now just export to environment (must be invoked w/ source, not bash)
export HIJETMETHODS=$internalHIJETMETHODS


FASTJETPATH=/home/cfmcginn/Packages/FastJet/fastjet-install
if [[ -d $FASTJETPATH ]]
then
    export FASTJETPATH=$FASTJETPATH
else
    echo "Given FASTJETPATH '$FASTJETPATH' does not exist. will not export to environment, some make may fail"   
fi
