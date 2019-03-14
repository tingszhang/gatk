#!/usr/bin/env bash

function getFileType()
{
  FT="UNSTRUCTURED_TEXT"
  bn=$( basename $1 )
  if [[ $bn =~ ^.*.bam$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="GENOME_READS"
  elif [[ $bn =~ ^.*.bam$ ]] ; then
    FT="EXOME_READS"
  elif [[ $bn =~ ^.*.bai$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="GENOME_READS_INDEX"
  elif [[ $bn =~ ^.*.bai$ ]] ; then
    FT="EXOME_READS_INDEX"
  elif [[ $bn =~ ^.*.vcf.gz$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf.gz$ ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf.gz.tbi$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF_INDEX"
  elif [[ $bn =~ ^.*.vcf.gz.tbi$ ]] ; then
    FT="SOMATIC_VCF_INDEX"
  elif [[ $bn =~ ^.*.vcf$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf$ ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf.idx$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF_INDEX"
  elif [[ $bn =~ ^.*.vcf.idx$ ]] ; then
    FT="SOMATIC_VCF_INDEX"
  fi

  echo "${FT}"
}

function getOutputFileID() 
{
  local fileBaseName=$1
  local runUUID=$2

  outputFileIndex=0
  if [[ "${runUUID}" == "84edd0ed-f084-47bb-af11-d8b47b9f1865" ]] ; then
    outputFileIndex=1
  elif [[ "${runUUID}" == "683dad15-3dea-4f35-8826-6d31f0e0c7bc" ]] ; then
    outputFileIndex=19
  elif [[ "${runUUID}" == "334ef0b9-ce97-49b3-8728-b6400396cde7" ]] ; then
    outputFileIndex=37
  elif [[ "${runUUID}" == "be5c2094-9ed2-4898-903d-cf519128ca48" ]] ; then
    outputFileIndex=55
  elif [[ "${runUUID}" == "5e7bc348-4745-4e6b-8231-bae0f57fc0b0" ]] ; then
    outputFileIndex=73
  elif [[ "${runUUID}" == "73a9ee75-6006-40f4-8f1b-681c38a501a8" ]] ; then
    outputFileIndex=91
  elif [[ "${runUUID}" == "662f5bfb-038c-491d-925b-896cc1038ff2" ]] ; then
    outputFileIndex=109
  elif [[ "${runUUID}" == "301fbc8e-2be1-41bd-847f-0dc4aff9f9af" ]] ; then
    outputFileIndex=127
  elif [[ "${runUUID}" == "001d8293-6c8c-41b9-8f07-274e49039d9b" ]] ; then
    outputFileIndex=145
  elif [[ "${runUUID}" == "1b557c70-a9e8-40f7-a4b2-d33c792c0255" ]] ; then
    outputFileIndex=163
  elif [[ "${runUUID}" == "88aba234-c27b-4bad-842e-9662704d64ca" ]] ; then
    outputFileIndex=181
  fi

  if [[ "${fileBaseName}" == "NexPond-359781" ]] ; then
    let outputFileIndex=$outputFileIndex+0
  elif [[ "${fileBaseName}" == "NexPond-359877" ]] ; then
    let outputFileIndex=$outputFileIndex+1
  elif [[ "${fileBaseName}" == "NexPond-360361" ]] ; then
    let outputFileIndex=$outputFileIndex+2
  elif [[ "${fileBaseName}" == "NexPond-360457" ]] ; then
    let outputFileIndex=$outputFileIndex+3
  elif [[ "${fileBaseName}" == "NexPond-361337" ]] ; then
    let outputFileIndex=$outputFileIndex+4
  elif [[ "${fileBaseName}" == "NexPond-361433" ]] ; then
    let outputFileIndex=$outputFileIndex+5
  elif [[ "${fileBaseName}" == "NexPond-362428" ]] ; then
    let outputFileIndex=$outputFileIndex+6
  elif [[ "${fileBaseName}" == "NexPond-363907" ]] ; then
    let outputFileIndex=$outputFileIndex+7
  elif [[ "${fileBaseName}" == "NexPond-445394" ]] ; then
    let outputFileIndex=$outputFileIndex+8
  fi

  echo $outputFileIndex
}

function createConcordanceJson() 
{
  echo "{"
  echo "  \"Concordance.gatk_docker\": \"broadinstitute/gatk-nightly:2019-02-26-4.1.0.0-31-g23bd0a2f8-SNAPSHOT\","
  echo ""
  echo "  \"Concordance.eval_vcf\": \"${1}\","
  echo "  \"Concordance.eval_vcf_idx\": \"${1}.idx\","
  echo "  \"Concordance.truth_vcf\": \"${2}\","
  echo "  \"Concordance.truth_vcf_idx\": \"${2}.tbi\","
  echo "  \"Concordance.intervals\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\","
  echo ""
  echo "  \"Concordance.default_disk_space_gb\": 512,"
  echo "  \"Concordance.mem_gb\": 32,"
  echo "  \"Concordance.boot_disk_size_gb\": 64"
  echo "}"
}

function createGenotypeConcordanceJson() 
{
  echo "{"
  echo "  \"GenotypeConcordanceTask.gatk_docker\": \"broadinstitute/gatk-nightly:2019-02-26-4.1.0.0-31-g23bd0a2f8-SNAPSHOT\","
  echo ""
  echo "  \"GenotypeConcordanceTask.call_vcf\": \"${1}\","
  echo "  \"GenotypeConcordanceTask.call_index\": \"${1}.idx\","
  echo "  \"GenotypeConcordanceTask.call_sample\": \"${1}.idx\","
  echo ""
  echo "  \"GenotypeConcordanceTask.truth_vcf\": \"${2}\","
  echo "  \"GenotypeConcordanceTask.truth_index\": \"${2}.tbi\","
  echo "  \"GenotypeConcordanceTask.truth_sample\": \"${2}.tbi\","
  echo ""
  echo "  \"GenotypeConcordanceTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\","
  echo ""
  echo "  \"GenotypeConcordanceTask.output_base_name\": \"${3}.tbi\","
  echo ""
  echo "  \"GenotypeConcordanceTask.default_disk_space_gb\": 512,"
  echo "  \"GenotypeConcordanceTask.mem_gb\": 32,"
  echo "  \"GenotypeConcordanceTask.boot_disk_size_gb\": 64"
  echo "}"
}

doToolDataCreation=false
doAnalysisDataCreation=true

if $doToolDataCreation ; then 

  outputFilesQueryFile=outputFilesQueryFile.mysql.txt
  testScenarioOutputFilesQueryFile=testScenarioOutputFilesQueryFile.mysql.txt

  echo "INSERT INTO DSPRegressionTesting.OutputFiles (fileType, path, timeCreated, md5sum, sourceType) VALUES " > ${outputFilesQueryFile}
  echo "INSERT INTO DSPRegressionTesting.TestScenarioOutputFiles (scenarioRun, outputFile) VALUES " > ${testScenarioOutputFilesQueryFile}

  scenarioRunIdBase=1
  scenarioId=1
  outputFileId=1
  isFirstFile=true
  for r in 84edd0ed-f084-47bb-af11-d8b47b9f1865 683dad15-3dea-4f35-8826-6d31f0e0c7bc 334ef0b9-ce97-49b3-8728-b6400396cde7 be5c2094-9ed2-4898-903d-cf519128ca48 5e7bc348-4745-4e6b-8231-bae0f57fc0b0 73a9ee75-6006-40f4-8f1b-681c38a501a8 662f5bfb-038c-491d-925b-896cc1038ff2 301fbc8e-2be1-41bd-847f-0dc4aff9f9af 001d8293-6c8c-41b9-8f07-274e49039d9b 1b557c70-a9e8-40f7-a4b2-d33c792c0255 88aba234-c27b-4bad-842e-9662704d64ca ; do

    echo "Processing run: ${r}"

    gatkDocker=$( cat ${r}/*SNAPSHOT*.json | jq '."ToolComparisonWdl.gatk_docker"'| tr -d '"' )

    tmpFile=$(mktemp)
    echo -n "  Getting outputs...  "
    ~/Development/cromshell/cromshell list-outputs $r 2>/dev/null | grep 'HaplotypeCaller' | grep -v 'timingInformation' | xargs gsutil ls -l | grep -v '^[ \t]*$'  | grep -v '^TOTAL' > $tmpFile
    echo 'DONE!'

    while read outFileInfo ; do

      echo "  Processing output file ${outputFileId}"

      outFile=$( echo $outFileInfo | awk '{print $3}' )

      shardIdx=$( echo $outFile | grep -o shard-[0-9]* | sed 's#shard-##g' )
      let scenarioRunId=$scenarioRunIdBase+$shardIdx

      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')

      FT="UNSTRUCTURED_TEXT"
      bn=$( basename $outFile )
      if [[ $bn =~ ^.*.bam$ ]] && [[ $bn =~ ^G.* ]] ; then
        FT="GENOME_READS"
      elif [[ $bn =~ ^.*.bam$ ]] ; then
        FT="EXOME_READS"
      elif [[ $bn =~ ^.*.bai$ ]] && [[ $bn =~ ^G.* ]] ; then
        FT="GENOME_READS_INDEX"
      elif [[ $bn =~ ^.*.bai$ ]] ; then
        FT="EXOME_READS_INDEX"
      elif [[ $bn =~ ^.*.vcf.gz$ ]] && [[ $bn =~ ^G.* ]] ; then
        FT="SOMATIC_VCF"
      elif [[ $bn =~ ^.*.vcf.gz$ ]] ; then
        FT="SOMATIC_VCF"
      elif [[ $bn =~ ^.*.vcf.gz.tbi$ ]] && [[ $bn =~ ^G.* ]] ; then
        FT="SOMATIC_VCF_INDEX"
      elif [[ $bn =~ ^.*.vcf.gz.tbi$ ]] ; then
        FT="SOMATIC_VCF_INDEX"
      elif [[ $bn =~ ^.*.vcf$ ]] && [[ $bn =~ ^G.* ]] ; then
        FT="SOMATIC_VCF"
      elif [[ $bn =~ ^.*.vcf$ ]] ; then
        FT="SOMATIC_VCF"
      elif [[ $bn =~ ^.*.vcf.idx$ ]] && [[ $bn =~ ^G.* ]] ; then
        FT="SOMATIC_VCF_INDEX"
      elif [[ $bn =~ ^.*.vcf.idx$ ]] ; then
        FT="SOMATIC_VCF_INDEX"
      fi

      md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )

      if ! $isFirstFile ; then
        echo -n ',' >> ${outputFilesQueryFile}
        echo -n ',' >> ${testScenarioOutputFilesQueryFile}
      fi

      echo "('$FT', '$outFile', '$timeStamp', '$md5sum', 'TOOL')" >> ${outputFilesQueryFile}
      echo "(${scenarioRunId}, ${outputFileId})" >> ${testScenarioOutputFilesQueryFile}
      let outputFileId=${outputFileId}+1
      isFirstFile=false
    done < $tmpFile

    let scenarioId=$scenarioId+1
    let scenarioRunIdBase=$scenarioRunIdBase+9

  done

  echo ";" >> ${outputFilesQueryFile}
  echo ';' >> ${testScenarioOutputFilesQueryFile}

  echo "Done processing TOOL output files."

fi

################################################################################

if $doAnalysisDataCreation ; then

  echo "Creating ANALYSIS file outputs."

  analysisRunsQueryFile=analysisRunsQueryFile.mysql.txt
  outputFilesAnalysisQueryFile=outputFilesAnalysisQueryFile.mysql.txt
  analysisOutputFilesQueryFile=analysisOutputFilesQueryFile.mysql.txt

  echo "INSERT INTO DSPRegressionTesting.AnalysisRuns (analysisInfoID, scenarioOutputForComparison, configuration) VALUES " > ${analysisRunsQueryFile}
  echo "INSERT INTO DSPRegressionTesting.OutputFiles (fileType, path, timeCreated, md5sum, sourceType) VALUES " > ${outputFilesAnalysisQueryFile}
  echo "INSERT INTO DSPRegressionTesting.AnalysisOutputFiles (analysis, outputFile) VALUES " > ${analysisOutputFilesQueryFile}

  outputFileId=199
  isFirstFile=true
  isFirstAnalysisFile=true
  analysisRunID=0
  for r in 84edd0ed-f084-47bb-af11-d8b47b9f1865 683dad15-3dea-4f35-8826-6d31f0e0c7bc 334ef0b9-ce97-49b3-8728-b6400396cde7 be5c2094-9ed2-4898-903d-cf519128ca48 5e7bc348-4745-4e6b-8231-bae0f57fc0b0 73a9ee75-6006-40f4-8f1b-681c38a501a8 662f5bfb-038c-491d-925b-896cc1038ff2 301fbc8e-2be1-41bd-847f-0dc4aff9f9af 001d8293-6c8c-41b9-8f07-274e49039d9b 1b557c70-a9e8-40f7-a4b2-d33c792c0255 88aba234-c27b-4bad-842e-9662704d64ca ; do

    echo "Processing run: ${r}"

    #################### 
    
    tmpFile=$(mktemp)
    echo -n "  Getting outputs..."
    ~/Development/cromshell/cromshell list-outputs $r 2>/dev/null | grep -v 'call-HaplotypeCallerTask/.*vcf[\.idx]*$' | xargs gsutil ls -l | grep -v '^[ \t]*$'  | grep -v '^TOTAL' > $tmpFile
    echo 'DONE!'

    echo -n "  Creating shard / input map..."
    shardInputMapFile=$(mktemp)
    grep 'HaplotypeCallerTask' $tmpFile | grep HaplotypeCallerTask | sed 's#.*shard-#shard-#g' | sed 's#\(.*NexPond-[0-9]*\).*#\1#g' | tr '/' '\t' > ${shardInputMapFile}
    echo 'DONE!'

    echo -n "  Sorting outputs into Concordance, GenotypeConcordance, and Timing sets"
    concordanceFile=$(mktemp)
    genotypeConcordanceFile=$(mktemp)
    timingFile=$(mktemp)
    grep 'call-Concordance' $tmpFile > $concordanceFile
    grep 'call-GenotypeConcordance' $tmpFile > $genotypeConcordanceFile
    grep 'call-HaplotypeCallerTask' $tmpFile > $timingFile 
    echo 'DONE!'

    echo "  Handling Concordance runs...  "
    lastShardIdx=999
    while read outFileInfo ; do
      echo "  Processing output file ${outputFileId}"
      
      outFile=$( echo $outFileInfo | awk '{print $3}' ) 
      shardIdx=$( echo $outFile | grep -o shard-[0-9]* | sed 's#shard-##g' )

      if [[ "${lastShardIdx}" != "${shardIdx}" ]] ; then
        let analysisRunID=$analysisRunID+1
        lastShardIdx=${shardIdx}
      fi

      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')
      truthFileBaseName=$(grep "shard-${shardIdx}" $shardInputMapFile | awk '{print $NF}' )
      inputFile="gs://broad-dsde-methods/cromwell-execution-36/ToolComparisonWdl/${r}/call-HaplotypeCallerTask/shard-${shardIdx}/${truthFileBaseName}.HC.vcf"
      truthFile="gs://broad-dsp-methods-regression-testing/inputData/${truthFileBaseName}.vcf.gz"
      originalOutFileID=$( getOutputFileID $truthFileBaseName $r )
      FT=$( getFileType $outFile )
      md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )
      
      if ! $isFirstFile ; then
        echo -n ',' >> ${analysisRunsQueryFile}
        echo -n ',' >> ${analysisOutputFilesQueryFile}
        echo -n ',' >> ${outputFilesAnalysisQueryFile}
      fi

      echo "(1, NULL, '$(createConcordanceJson $inputFile $truthFile | jq -c .)')" >> $analysisRunsQueryFile 

      grep -q "$analysisRunID, $originalOutFileID" $analysisOutputFilesQueryFile
      [ $? -ne 0 ] && echo "($analysisRunID, $originalOutFileID)" >> $analysisOutputFilesQueryFile
      
      echo "('$FT', '$outFile', '$timeStamp', '$md5sum', 'ANALYSIS')" >> ${outputFilesAnalysisQueryFile}
      echo "($analysisRunID, $outputFileId)" >> $analysisOutputFilesQueryFile

      let outputFileId=${outputFileId}+1
      isFirstFile=false
    done < $concordanceFile

    exit

    #################### 
    echo "Handling Concordance files..."
    while read outFileInfo; do 
      echo "  Processing output file ${outputFileId}"
      outFile=$( echo $outFileInfo | awk '{print $3}' ) 
      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')
      FT=$( getFileType $outFile )
      md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )
      if ! $isFirstFile ; then
        echo -n ',' >> ${outputFilesAnalysisQueryFile}
        echo -n ',' >> ${analysisOutputFilesQueryFile}
      fi
      echo "('$FT', '$outFile', '$timeStamp', '$md5sum', 'ANALYSIS')" >> ${outputFilesAnalysisQueryFile}
      echo "(NULL, ${outputFileId})" >> ${analysisOutputFilesQueryFile}
      let outputFileId=${outputFileId}+1
      isFirstFile=false
    done < $concordanceFile 

    echo "Handling GenotypeConcordance files..."
    while read outFileInfo; do 
      echo "  Processing output file ${outputFileId}"
      outFile=$( echo $outFileInfo | awk '{print $3}' ) 
      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')
      FT=$( getFileType $outFile )
      md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )
      echo ",('$FT', '$outFile', '$timeStamp', '$md5sum', 'ANALYSIS')" >> ${outputFilesAnalysisQueryFile}

      echo ",(NULL, ${outputFileId})" >> ${analysisOutputFilesQueryFile}

      let outputFileId=${outputFileId}+1
    done < $genotypeConcordanceFile

    echo "Handling Timing files..."
    while read outFileInfo; do 
      echo "  Processing output file ${outputFileId}"
      outFile=$( echo $outFileInfo | awk '{print $3}' ) 
      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')
      FT=$( getFileType $outFile )
      md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )
      echo ",('$FT', '$outFile', '$timeStamp', '$md5sum', 'ANALYSIS')" >> ${outputFilesAnalysisQueryFile}

      echo ",(NULL, ${outputFileId})" >> ${analysisOutputFilesQueryFile}

      let outputFileId=${outputFileId}+1
    done < $timingFile

  done

  echo ";" >> ${outputFilesAnalysisQueryFile}
  echo ';' >> ${analysisOutputFilesQueryFile}

  echo "Done processing ANALYSIS output files."

fi


