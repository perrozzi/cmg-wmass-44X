#!/bin/bash

# use_PForNoPUorTKmet=( 2 1 0 ) # 0:PF, 1:NOPU, 2:TK 
use_PForNoPUorTKmet=( 2 ) # 0:PF, 1:NOPU, 2:TK 

# MOMENTUM CORRECTION VARIATIONS
momcorr_scale_variations=( 0 )

# RECOIL CORRECTION VARIATIONS
# Recoil_resol_variations=( 0 1 -1 0  0 )
# Recoil_scale_variations=( 0 0  0 1 -1 )
Recoil_resol_variations=( 0 )
Recoil_scale_variations=( 0 )
# Recoil_resol_variations=( 0 1 -1 )
# Recoil_scale_variations=( 0 0  0 )
# Recoil_resol_variations=( 0 0 )
# Recoil_scale_variations=( 0 1 )
useRecoilCorr=( 1 ) # 1: YES, 0: NO
# useRecoilCorr=( 1 ) # 1: YES, 0: NO

# infile_run='launch_analysis.py'
infile_run='launch_analysis_testSplit_bash.py'

# echo ${#Recoil_resol_variations[@]}

run_all_or_just_fit=3 #  2 = W/Z ANALYSIS, 1 = RUN ALL,  0 = RUN FIT ONLY, 3 = W/Z ANALYSIS + QCD FIT
run_W_or_Z=0 #  0 = W,  1 = Z,  2 = both (W and Z)
fit_W_or_Z="W,Z" # "W" or "Z" or "W,Z"

for ((j=0; j < ${#use_PForNoPUorTKmet[@]} ; j++));
    do
    
  for ((k=0; k < ${#useRecoilCorr[@]} ; k++));
      do
      
    for ((i=0; i < ${#Recoil_resol_variations[@]} ; i++));
        do
        
      for ((m=0; m < ${#Recoil_scale_variations[@]} ; m++));
          do
          
        for ((h=0; h < ${#momcorr_scale_variations[@]} ; h++));
            do
            
            if [ ${useRecoilCorr[${k}]} -eq 0 ] ; then
              if [ ${Recoil_resol_variations[${i}]} -ne 0 ] || [ ${Recoil_scale_variations[${m}]} -ne 0 ]; then
                continue
              fi
            fi
            
            echo "use_PForNoPUorTKmet " ${use_PForNoPUorTKmet[${j}]} " (j="$j")"
            echo "useRecoilCorr " ${useRecoilCorr[${k}]} " (k="$k")"
            echo "Recoil_resol_variations " ${Recoil_resol_variations[${i}]} " (i="$i")"
            echo "Recoil_scale_variations " ${Recoil_scale_variations[${m}]} " (m="$m")"
            echo "momcorr_scale_variations " ${momcorr_scale_variations[${h}]} " (h="$h")"

            # echo ${Recoil_resol_variations[${i}]} ${Recoil_scale_variations[${m}]}; 
            # for (resol,scale) in Recoil_scale_variations:
            # echo "i="$i " Recoil_resol_variations " ${Recoil_resol_variations[${i}]} " Recoil_scale_variations "${Recoil_scale_variations[${m}]}
            # echo "h="$h " momcorr_scale_variations " ${momcorr_scale_variations[${h}]}
            cp launch_analysis_testSplit_bkp_bash.py $infile_run

            sed -i "s/.*fit_W_or_Z = \"/fit_W_or_Z = \"$fit_W_or_Z\" \# \"W\" or \"Z\"/g" $infile_run
            sed -i "s/.*use_PForNoPUorTKmet = .*/use_PForNoPUorTKmet = ${use_PForNoPUorTKmet[${j}]}\ \#\ 0\:PF\,\ 1\:NOPU\,\ 2\:TK\;/g" $infile_run
            sed -i "s/.*useRecoilCorr = .*/useRecoilCorr = ${useRecoilCorr[${k}]}\ \#\ 1\:YES\,\ 0\:NO/g" $infile_run
            
            echo " "
            
            ###################################
            echo '## RUN W AND Z ANALYSES ##'
            # ###################################

            sed -i "s/.*RecoilCorrResolutionNSigma =.*/RecoilCorrResolutionNSigma = \"${Recoil_resol_variations[${i}]}\";/" $infile_run
            sed -i "s/.*RecoilCorrScaleNSigma =.*/RecoilCorrScaleNSigma = \"${Recoil_scale_variations[${m}]}\";/" $infile_run
            sed -i "s/.*GlobalSmearingRochCorrNsigma =.*/GlobalSmearingRochCorrNsigma = ${momcorr_scale_variations[${h}]};/" $infile_run
            grep RecoilCorrResolutionNSigma\ = $infile_run
            grep RecoilCorrScaleNSigma\ = $infile_run
            grep GlobalSmearingRochCorrNsigma\ = $infile_run
            
            # continue
            
            if [ $run_W_or_Z -eq 0 ]; then # run W
              sed -i 's/runWanalysis = 0;/runWanalysis = 1;/g' $infile_run
              sed -i 's/runZanalysis = 1;/runZanalysis = 0;/g' $infile_run
            elif [ $run_W_or_Z -eq 1 ]; then # run Z
              sed -i 's/runWanalysis = 1;/runWanalysis = 0;/g' $infile_run
              sed -i 's/runZanalysis = 0;/runZanalysis = 1;/g' $infile_run
            else
              sed -i 's/runWanalysis = 0;/runWanalysis = 1;/g' $infile_run
              sed -i 's/runZanalysis = 0;/runZanalysis = 1;/g' $infile_run
            fi
            sed -i 's/mergeSigEWKbkg = 1;/mergeSigEWKbkg = 0;/g' $infile_run
            sed -i 's/runR_WdivZ= 1;/runR_WdivZ= 0;/g' $infile_run
            sed -i 's/run_BuildSimpleTemplates= 1;/run_BuildSimpleTemplates= 0;/g' $infile_run
            sed -i 's/runPrepareDataCards = 1;/runPrepareDataCards = 0;/g' $infile_run
            sed -i 's/runClosureTestLikeLihoodRatioAnsMergeResults = 1;/runClosureTestLikeLihoodRatioAnsMergeResults = 0;/g' $infile_run
            sed -i 's/runWSigBkgFit = 1;/runWSigBkgFit = 0;/g' $infile_run
            if [ $run_all_or_just_fit -eq 3 ]; then
              sed -i 's/controlplots = 0;/controlplots = 1;/g' $infile_run
            fi
            
            if [ $run_all_or_just_fit -gt 0 ]; then
              # echo "not launching"
              python $infile_run
            fi
            
            ##############################################################
            echo '## CHECK IF RUN W AND Z ANALYSES ARE STILL RUNNING ##'
            ##############################################################
            still_running=1
            while [ $still_running -eq 1 ]; do
              line=$(ps aux | grep perrozzi | grep analysis |grep root)
              # echo $line
              if [[ ! ${line} =~ "runWanalysis" ]]; then
                if [[ ! ${line} =~ "runZanalysis" ]]; then
                  still_running=0
              # if still_running==1: print 'runWanalysis or runZanalysis still running'
              # else: print 'runWanalysis and runZanalysis not running anymore'
                fi
              else
              # else: 
                sleep 10
                # print 'runWanalysis or runZanalysis still running',line
              fi
            done

            ###################################
            echo '## MERGE EWK AND TT     ##'
            ###################################
            sed -i 's/runWanalysis = 1;/runWanalysis = 0;/g' $infile_run
            sed -i 's/runZanalysis = 1;/runZanalysis = 0;/g' $infile_run
            sed -i 's/mergeSigEWKbkg = 0;/mergeSigEWKbkg = 1;/g' $infile_run
            sed -i 's/runR_WdivZ= 1;/runR_WdivZ= 0;/g' $infile_run
            sed -i 's/run_BuildSimpleTemplates= 1;/run_BuildSimpleTemplates= 0;/g' $infile_run
            sed -i 's/runPrepareDataCards = 1;/runPrepareDataCards = 0;/g' $infile_run
            sed -i 's/runClosureTestLikeLihoodRatioAnsMergeResults = 1;/runClosureTestLikeLihoodRatioAnsMergeResults = 0;/g' $infile_run
            sed -i 's/runWSigBkgFit = 1;/runWSigBkgFit = 0;/g' $infile_run
            if [ $run_all_or_just_fit -gt 0 ]; then
              # echo "not launching"
              python $infile_run
            fi

              
            ###################################
            echo '## FIT QCD BACKGROUND   ##'
            ###################################
            sed -i 's/runWanalysis = 1;/runWanalysis = 0;/g' $infile_run
            sed -i 's/runZanalysis = 1;/runZanalysis = 0;/g' $infile_run
            sed -i 's/mergeSigEWKbkg = 1;/mergeSigEWKbkg = 0;/g' $infile_run
            sed -i 's/runR_WdivZ= 1;/runR_WdivZ= 0;/g' $infile_run
            sed -i 's/run_BuildSimpleTemplates= 1;/run_BuildSimpleTemplates= 0;/g' $infile_run
            sed -i 's/runPrepareDataCards = 1;/runPrepareDataCards = 0;/g' $infile_run
            sed -i 's/runClosureTestLikeLihoodRatioAnsMergeResults = 1;/runClosureTestLikeLihoodRatioAnsMergeResults = 0;/g' $infile_run
            sed -i 's/runWSigBkgFit = 0;/runWSigBkgFit = 1;/g' $infile_run
            sed -i 's/controlplots = 0;/controlplots = 1;/g' $infile_run
            if [ $run_all_or_just_fit -eq 3 ]; then
              # echo "not launching"
              python $infile_run
            fi
              
            ###################################
            echo '## RUN R = W div Z      ##'
            ###################################
            sed -i 's/runWanalysis = 1;/runWanalysis = 0;/g' $infile_run
            sed -i 's/runZanalysis = 1;/runZanalysis = 0;/g' $infile_run
            sed -i 's/mergeSigEWKbkg = 1;/mergeSigEWKbkg = 0;/g' $infile_run
            sed -i 's/runR_WdivZ= 0;/runR_WdivZ= 1;/g' $infile_run
            sed -i 's/run_BuildSimpleTemplates= 1;/run_BuildSimpleTemplates= 0;/g' $infile_run
            sed -i 's/runPrepareDataCards = 1;/runPrepareDataCards = 0;/g' $infile_run
            sed -i 's/runClosureTestLikeLihoodRatioAnsMergeResults = 1;/runClosureTestLikeLihoodRatioAnsMergeResults = 0;/g' $infile_run
            sed -i 's/runWSigBkgFit = 1;/runWSigBkgFit = 0;/g' $infile_run
            # python $infile_run
              
            #########################################
            echo '## RUN BUILD SIMPLE TEMPLATES ##'
            #########################################
            sed -i 's/runWanalysis = 1;/runWanalysis = 0;/g' $infile_run
            sed -i 's/runZanalysis = 1;/runZanalysis = 0;/g' $infile_run
            sed -i 's/mergeSigEWKbkg = 1;/mergeSigEWKbkg = 0;/g' $infile_run
            sed -i 's/runR_WdivZ= 1;/runR_WdivZ= 0;/g' $infile_run
            sed -i 's/run_BuildSimpleTemplates= 0;/run_BuildSimpleTemplates= 1;/g' $infile_run
            sed -i 's/runPrepareDataCards = 1;/runPrepareDataCards = 0;/g' $infile_run
            sed -i 's/runClosureTestLikeLihoodRatioAnsMergeResults = 1;/runClosureTestLikeLihoodRatioAnsMergeResults = 0;/g' $infile_run
            sed -i 's/runWSigBkgFit = 1;/runWSigBkgFit = 0;/g' $infile_run
            # python $infile_run

            ####################################
            echo '## RUN PREPARE DATACARDS ##'
            ####################################
            sed -i 's/runWanalysis = 1;/runWanalysis = 0;/g' $infile_run
            sed -i 's/runZanalysis = 1;/runZanalysis = 0;/g' $infile_run
            sed -i 's/mergeSigEWKbkg = 1;/mergeSigEWKbkg = 0;/g' $infile_run
            sed -i 's/runR_WdivZ= 1;/runR_WdivZ= 0;/g' $infile_run
            sed -i 's/run_BuildSimpleTemplates= 1;/run_BuildSimpleTemplates= 0;/g' $infile_run
            sed -i 's/runPrepareDataCards = 0;/runPrepareDataCards = 1;/g' $infile_run
            sed -i 's/runClosureTestLikeLihoodRatioAnsMergeResults = 1;/runClosureTestLikeLihoodRatioAnsMergeResults = 0;/g' $infile_run
            sed -i 's/runWSigBkgFit = 1;/runWSigBkgFit = 0;/g' $infile_run
            if [ $run_all_or_just_fit -lt 2 ]; then
              # echo "not launching"
              python $infile_run
            fi

            ###################################################
            echo '## RUN LIKELIHOOD FIT AND MERGE RESULTS ##'
            ###################################################
            sed -i 's/runWanalysis = 1;/runWanalysis = 0;/g' $infile_run
            sed -i 's/runZanalysis = 1;/runZanalysis = 0;/g' $infile_run
            sed -i 's/mergeSigEWKbkg = 1;/mergeSigEWKbkg = 0;/g' $infile_run
            sed -i 's/runR_WdivZ= 1;/runR_WdivZ= 0;/g' $infile_run
            sed -i 's/run_BuildSimpleTemplates= 1;/run_BuildSimpleTemplates= 0;/g' $infile_run
            sed -i 's/runPrepareDataCards = 1;/runPrepareDataCards = 0;/g' $infile_run
            
            sed -i 's/runClosureTestLikeLihoodRatioAnsMergeResults = 0;/runClosureTestLikeLihoodRatioAnsMergeResults = 1;/g' $infile_run
            sed -i 's/mergeResults = 1;/mergeResults = 0;/g' $infile_run
            sed -i 's/runWSigBkgFit = 1;/runWSigBkgFit = 0;/g' $infile_run
            # sed -i 's/runClosureTestLikeLihoodRatioAnsMergeResults = 1;/runClosureTestLikeLihoodRatioAnsMergeResults = 0;/g' $infile_run
            # sed -i 's/mergeResults = 0;/mergeResults = 1;/g' $infile_run
            
            if [ $run_all_or_just_fit -lt 2 ]; then
              # echo "not launching"
              python $infile_run
            fi

            echo " "
        done
      done
    done
  done
done