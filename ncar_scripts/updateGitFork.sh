#!/usr/bin/env bash
##################################################################
# This script is intended to update our release version of
# clubb (the one "publicly" visible). It should be ran nightly.
#
# Input options:
#    * source repo: URL for the original repo (the repo yours was forked from)
#    * destination repo: URL for your repo
#    * source branch: branch you want to pull updates from (defaults to master)
#    * destination branch: branch you want to push updates to (defaults to master)
#    * force update: adds the '-f' flag to the push command, thus force pushes
#    * mirror (-m): update all branches (as long as they have the same name)
#
# More usage info can be found here:
# https://github.com/larson-group/sys_admin/wiki/GitUpdateScripts#updategitforksh
#
# Author: Nicolas Strike, 2018
##################################################################
workdir=/tmp/temprepo
srcbranch=main
destbranch=main
mirror=""
forceupdate=""
email=0
allbranches=0
tags=""
maxNetworkAttempts=15
internetConnectivityTestUrl="github.com"
    while [ "$1" != "" ]; do
        case $1 in
            -s | --source-repo )    shift
                                    srcrepo=$1
                                    ;;
            -d | --destination-repo )     shift
                                    destrepo=$1
                                    ;;
            --source-branch )       shift
                                    srcbranch=$1
                                    ;;
            --dest-branch )         shift
                                    destbranch=$1
                                    ;;
            -f | --force-update )   forceupdate="-f"
                                    ;;
            -m | --mirror )         destbranch="--mirror"
                                    ;;
            -a | --all-branches )   allbranches=1
                                    ;;
            -e | --enable-email )   email=1
                                    ;;
            -t | --tags )           tags="--tags"
                                    ;;
            -h | --help )           usage
                                    exit
                                    ;;
        esac
        shift
done


networkSafeOperation(){
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

rm -rf $workdir
networkSafeOperation git clone $srcrepo $workdir --mirror

cd $workdir
git remote add destination_repo $destrepo

networkSafeOperation git fetch $tags
echo "Pushing updates from $srcrepo:$srcbranch to $destrepo:$destbranch"
if [ $allbranches != 0 ]
then
    # Alert user a parameter will be ignored #
    if [ $destbranch == "--mirror" ]
    then
        echo "The mirror is not compatible with the allbranches parameter, it will be ignored. All branches will still be pushed using a for-each branch loop instead of using git mirror. The difference is branches on the destination repo will not be deleted, whereas the would've been if --mirror was used."
    fi
    if [ $srcbranch != "master" ]
    then
        echo "The source branch parameter is being ignored because --all-branches was passed. All branches on the source repo will be pushed to the same branch on the destination repo."
    fi
    # ----------------------------------------
    # Finally push the updates
    for branch in $(git for-each-ref --format='%(refname)' refs/heads/); do
        networkSafeOperation git push $forceupdate -u destination_repo $branch $tags
    done
else
    networkSafeOperation git push $forceupdate -u destination_repo refs/heads/$srcbranch:refs/heads/$destbranch $tags
fi
errorsDetected=$?
wait


cd -
rm -rf $workdir

#####################################################################################################
#This sends out an email notification to system admins if any errors were detected during runtime.
# The content of the email can be found in $emailMessage
#####################################################################################################
sendEmailNotification(){

	if [ $errorsDetected != 0 ]
	then
	  echo "Errors during run detected. Sending email notification."
	  #Set the email
	  EMAIL=cmille73@ucar.edu
          echo -e "Hello,\nWhen attempting to automatically update $destrepo:$destbranch from $srcrepo:$srcbranch there was an issue. Please run the updateGitFork.sh script manually. Additional troubleshooting info can be found here: https://github.com/larson-group/sys_admin/wiki/GitUpdateScripts \n\n`date +%Y-%m-%d` fork autoupdate script `hostname`\n\nThank you,\n  Team Nightly" | mail -s "Fork autoupdate script errors: `hostname`" $EMAIL
	  exit 1
	else
	   echo "No errors detected. Run successful."
	fi
}
if [ $email -eq 1 ]
then
	sendEmailNotification
fi
