# tTEscanR Function Map <a href="https://your-link-here"><img src="https://github.com/user-attachments/assets/8cca530a-f0a5-4284-bf2e-cf030a8193fa" alt="logo_tTEscanR" align="right" width="120"></a>
<img src="https://img.shields.io/badge/Language-R-blue.svg" style="zoom:100%;" />

<table>
  <tr>
    <th>Section</th>
    <th>File</th>
    <th>Function</th>
    <th>Function Status</th>
  </tr>

  <!-- tTEscanRObject -->
  <tr>
    <td rowspan="5">tTEscanR Object</td>
    <td rowspan="5">tTEscanRObject.R</td>
    <td>createObject( )</td>
    <td rowspan="2">Visible</td>
  </tr>
  <tr>
    <td>updateObject( )</td>
  </tr>
  <tr>
    <td>defineNewData( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>checkInputCombinations( )</td>
  </tr>
  <tr>
    <td>isInObject( )</td>
  </tr>

  <!-- General pipeline -->
  <tr>
    <td>General Pipeline</td>
    <td>tTEscanR_General.R</td>
    <td>runPipeline( )</td>
    <td>Visible</td>
  </tr>

  <!-- Codon Usage -->
  <tr>
    <td rowspan="8">Codon Usage</td>
    <td rowspan="2">CodonUsage.R</td>
    <td>computeCodonUsage( )</td>
    <td rowspan="2">Visible</td>
  </tr>
  <tr>
    <td>computeAnticodonUsage( )</td>
  </tr>
  <tr>
    <td rowspan="6">tTE_HelperFunctions_CodonUsage.R</td>
    <td>computeCodonExonicBackground( )</td>
    <td rowspan="4">Visible</td>
  </tr>
  <tr>
    <td>computeCorrelationBackground( )</td>
  </tr>
  <tr>
    <td>computeMeanUsage( )</td>
  </tr>
  <tr>
    <td>computeMetricsCodonUsage( )</td>
  </tr>
  <tr>
    <td>consistencyWithCodonFreq( )</td>
    <td rowspan="2">Internal</td>
  </tr>
  <tr>
    <td>CheckNames_tRNAs( )</td>
  </tr>

  <!-- Codon Freq Per Gene -->
  <tr>
    <td rowspan="7">Codon Freq Per Gene</td>
    <td>CodonFreqPerGene.R</td>
    <td>getCodonFreq( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td rowspan="6">tTE_HelperFunctions_CodonFreq.R</td>
    <td>extractCodons( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>extractGenes( )</td>
    <td rowspan="5">Internal</td>
  </tr>
  <tr>
    <td>extractSequences( )</td>
  </tr>
  <tr>
    <td>callingEnsembl( )</td>
  </tr>
  <tr>
    <td>checkFASTAFormat( )</td>
  </tr>
  <tr>
    <td>fromFASTAtoTable( )</td>
  </tr>

  <!-- Codon Pool Contribution -->
  <tr>
    <td rowspan="4">Pool Contribution</td>
    <td rowspan="4">ExamineCodonPoolContribution.R</td>
    <td>showPoolContribution( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>computeIndividualGeneCorrelation( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>analyzeTopGeneImpact( )</td>
  </tr>
  <tr>
    <td>getOrCompute( )</td>
  </tr>

  <!-- Amino Acid Usage -->
  <tr>
    <td rowspan="3">Amino Acid Usage</td>
    <td rowspan="3">AAUsage.R</td>
    <td>computeAAUsage( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>groupAA( )</td>
    <td rowspan="2">Internal</td>
  </tr>
  <tr>
    <td>retrieveAAUsageData( )</td>
  </tr>

  <!-- Theoretical Translation Efficiency -->
  <tr>
    <td rowspan="4">Theoretical Translation Efficiency</td>
    <td rowspan="4">TheoreticalTranslationEfficiency.R</td>
    <td>computeTheoreticalTE( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>computeStatisticalSignificance( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>computeCorrelation( )</td>
  </tr>
  <tr>
    <td>filterMatrix( )</td>
  </tr>

  <!-- Abundance analysis -->
  <tr>
    <td rowspan="13">Differential Expression Analysis</td>
    <td rowspan="3">DifferentialExpressionAnalysis.R</td>
    <td>runDEAnalysis( )</td>
    <td rowspan="3">Visible</td>
  </tr>
  <tr>
    <td>computeDEResults( )</td>
  </tr>
  <tr>
    <td>plotDEResults( )</td>
  </tr>
  <tr>
    <td rowspan="10">tTE_HelperFunctions_DEA.R</td>
    <td>computeSizeCorrection( )</td>
    <td rowspan="10">Internal</td>
  </tr>
  <tr>
    <td>computeAllPairwiseComp( )</td>
  </tr>
  <tr>
    <td>computeDESeq2( )</td>
  </tr>
  <tr>
    <td>targetedApproach( )</td>
  </tr>
  <tr>
    <td>makeScatterPlot( )</td>
  </tr>
  <tr>
    <td>runDimReduct( )</td>
  </tr>
  <tr>
    <td>generateVolcanoPlot( )</td>
  </tr>
  <tr>
    <td>heatmapFontSize( )</td>
  </tr>
  <tr>
    <td>produceHeatmapDiffExp( )</td>
  </tr>
  <tr>
    <td>produceElbowPlot( )</td>
  </tr>

  <!-- Permutation Test -->
  <tr>
    <td rowspan="2">Permutation Test</td>
    <td rowspan="2">PermutationTest.R</td>
    <td>getPermutationDist( )</td>
    <td rowspan="2">Visible</td>
  </tr>
  <tr>
    <td>obtainSignificance( )</td>
  </tr>

  <!-- General Functions -->
  <tr>
    <td rowspan="26">General Functions</td>
    <td rowspan="7">GeneralProcessing.R</td>
    <td>tRNAGetMatrix( )</td>
    <td rowspan="7">Visible</td>
  </tr>
  <tr>
    <td>tRNASetCutoff( )</td>
  </tr>
  <tr>
    <td>tRNAFilterCuts( )</td>
  </tr>
  <tr>
    <td>tRNASetGenes( )</td>
  </tr>
  <tr>
    <td>transformFormat( )</td>
  </tr>
  <tr>
    <td>mergeMatrices( )</td>
  </tr>
  <tr>
    <td>groupConditions( )</td>
  </tr>
  <tr>
    <td rowspan="9">tTE_HelperFunctions_GeneralProcessing.R</td>
    <td>transformCounts( )</td>
    <td rowspan="9">Internal</td>
  </tr>
  <tr>
    <td>cutoffMatrix( )</td>
  </tr>
  <tr>
    <td>filteringCutoffs( )</td>
  </tr>
  <tr>
    <td>computeCorrelations( )</td>
  </tr>
  <tr>
    <td>selectionCutoff( )</td>
  </tr>
  <tr>
    <td>iteratetRNACutoff( )</td>
  </tr>
  <tr>
    <td>correlationCutoffPlot( )</td>
  </tr>
  <tr>
    <td>selectionCutoffPlot( )</td>
  </tr>
  <tr>
    <td>cumulativeSelectionCutoffPlot( )</td>
  </tr>
  <tr>
    <td rowspan="6">tTE_HelperFunctions_GeneralInternal.R</td>
    <td>selectDefaultData( )</td>
    <td rowspan="6">Internal</td>
  </tr>
  <tr>
    <td>identifyInputFormat( )</td>
  </tr>
  <tr>
    <td>filterByMetadata( )</td>
  </tr>
  <tr>
    <td>checkCodonFreqTable( )</td>
  </tr>
  <tr>
    <td>checkDataFrame( )</td>
  </tr>
  <tr>
    <td>checkGeneAnnotation( )</td>
  </tr>
  <tr>
    <td rowspan="4">tTE_HelperFunctions_GeneralTranslation.R</td>
    <td>featuresToAA( )</td>
    <td>Visible</td>
  </tr>
  </tr>
    <td>performTranslation( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>retrieveTranslation( )</td>
  </tr>
  <tr>
    <td>isEnsemblID( )</td>
  </tr>

  <!-- Plots -->
<tr>
    <td rowspan="22">Plots</td>
    <td rowspan="7">tTEscanPlots.R</td>
    <td>plotDistribution( )</td>
    <td rowspan="6">Visible</td>
</tr>
<tr>
    <td>plotTargetComparison( )</td>
</tr>
<tr>
    <td>plotProportion( )</td>
</tr>
<tr>
    <td>plotScores( )</td>
</tr>
<tr>
    <td>plotCorrelation( )</td>
</tr>
<tr>
    <td>plotPermutation( )</td>
</tr>
<tr>
    <td>correlationCutoffPlot( )</td>
    <td>Internal</td>
  </tr>

<tr>
    <td rowspan="15">tTE_HelperFunctions_tTEscanRPlots.R</td>
    <td>savePlot( )</td>
    <td rowspan="15">Internal</td>
</tr>
<tr>
    <td>getSafeColorScale( )</td>
</tr>
<tr>
    <td>getAnnotData( )</td>
</tr>
<tr>
    <td>checkValueInData( )</td>
</tr>
<tr>
    <td>checkDataInLongFormat( )</td>
</tr>
<tr>
    <td>getOutputName( )</td>
</tr>
<tr>
    <td>generateDistPlot( )</td>
</tr>
<tr>
    <td>drawBarCountsPlot( )</td>
</tr>
<tr>
    <td>drawDonutPlot( )</td>
</tr>
<tr>
    <td>computeRings( )</td>
</tr>
<tr>
    <td>drawRadarPlot( )</td>
</tr>
<tr>
    <td>generateProportionPlot( )</td>
</tr>
<tr>
    <td>significanceSymbol( )</td>
</tr>
<tr>
    <td>computeTEsignificance( )</td>
</tr>
<tr>
    <td>computeBoxplotSignificance( )</td>
</tr>

</table>
