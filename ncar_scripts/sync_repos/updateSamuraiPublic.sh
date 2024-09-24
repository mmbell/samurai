#!/usr/bin/env bash

##################################################################
# This script is intended to update our release version of
# samurai (the one "publicly" visible). It should be ran nightly.
#
# Author: Nicolas Strike, 2019
##################################################################

# To add a new branch from samurai-dev to samurai (release), simply add it to the array here and it will be added overnight and updated automatically. Note that there must be a space between items and the parethesis.
branches=( main )

# List of patterns describing lfs files to untrack in git-lfs. This prevents them from being automatically downloaded when a client clones the repo, but leaves them available to download if desired (with `git lfs pull`)
undesiredLfsPatterns=( )


maxNetworkAttempts=15
internetConnectivityTestUrl="github.com"

# Need to fill in name of private development samurai repo (Probably git@github.com:NCAR/samurai-dev.git)

# UNCOMMENT THE FOLLOWING LINE AND VERIFY THE REPO NAMES BEFORE RUNNING
#updateSrc=git@github.com:NCAR/samurai-dev.git
updateDest=git@github.com:mmbell/samurai.git

updateBranch() {
    bash updateGitFork.sh --source-repo $updateSrc --destination-repo $updateDest --source-branch $1 --dest-branch $1 -f -t
}

networkSafeOperation (){
  #limited to 8 arguments, see below

  tryAgain="true"
  netAttempts=$maxNetworkAttempts
  # initial wait time before the second network attempt in seconds.
  # Sequential attempts get a x2 multiplier
  netWaitDelay=2
  while [ $tryAgain == "true" -a $netAttempts -gt 0 ]; do

    #Run the given command with up to 8 arguments
    $1 $2 $3 $4 $5 $6 $7 $8 $9

    cloneCode=$?
    netAttempts=$((netAttempts-1))
    #check if the repo successfully cloned
    if [ $cloneCode == 0 ]
    then
        tryAgain="false"
    else
        echo "Testing network connectivity."
        ping -c 4 $internetConnectivityTestUrl
        pingCode=$?
        if [ $pingCode != 0 ]
        then
          echo "Failed to perform network operation '$1 $2 $3 $4 $5 $6 $7' due to a network issue. $netAttempts attempts left, waiting $netWaitDelay seconds."
          sleep $netWaitDelay
          netWaitDelay=$((netWaitDelay*2))
        else
          echo "Failed to perform operation '$1 $2 $3 $4 $5 $6 $7', but no network issue was detected. This must be handeled manually and is likely a git conflict."
          tryAgain="false"
        #pingCode
        fi
    #cloneCode
    fi
    done
}

updateVersionFiles(){
    rm -rf /tmp/temp_release_repo
    networkSafeOperation git clone $updateDest /tmp/temp_release_repo
    cd /tmp/temp_release_repo
        git log -10 > src/CLUBB_core/version_clubb_core.txt
        git log -10 > src/SILHS/version_silhs.txt
        git add src/CLUBB_core/version_clubb_core.txt src/SILHS/version_silhs.txt
        git commit -a -m "Updating the clubb core and silhs version files"
        networkSafeOperation git push
    cd -
    rm -rf /tmp/temp_release_repo
}

untrackLfs(){
    for i in "${undesiredLfsPatterns[@]}"
    do
	    git lfs untrack $i
    done
}

updateAllBranches(){
    for i in "${branches[@]}"
    do
	    updateBranch $i
    done
}

printf "\n\n######################################\nUpdating repo at $updateDest from $updateSrc\n######################################\n"

    printf "\n###### Updating all branches ######\n"
        updateAllBranches
    printf "Done updating branches.\n"

    #printf "\n###### Updating untracking lfs files ######\n"
    #    untrackLfs
    #printf "Done untracking lfs files.\n"

    #printf "\n###### Updating version files ######\n"
    #    updateVersionFiles
    #printf "Done updating version files\n"

printf "\n###### DONE ######\n"
