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
    <td>Create_tTEscanR_Object( )</td>
    <td rowspan="2">Visible</td>
  </tr>
  <tr>
    <td>Update_tTEscanR_Object( )</td>
  </tr>
  <tr>
    <td>DefineNewData( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>CheckInputCombinations( )</td>
  </tr>
  <tr>
    <td>IsIn_tTEscanR_Object( )</td>
  </tr>

  <!-- General pipeline -->
  <tr>
    <td>General Pipeline</td>
    <td>tTEscanR_General.R</td>
    <td>Run_tTEscanR_pipeline( )</td>
    <td>Visible</td>
  </tr>

  <!-- Codon Usage -->
  <tr>
    <td rowspan="8">Codon Usage</td>
    <td rowspan="2">CodonUsage.R</td>
    <td>ComputeCodonUsage( )</td>
    <td rowspan="2">Visible</td>
  </tr>
  <tr>
    <td>ComputeAnticodonUsage( )</td>
  </tr>
  <tr>
    <td rowspan="6">tTE_HelperFunctions_CodonUsage.R</td>
    <td>ComputeCodonExonicBackground( )</td>
    <td rowspan="4">Visible</td>
  </tr>
  <tr>
    <td>ComputeCorrelationBackground( )</td>
  </tr>
  <tr>
    <td>ComputeMeanUsage( )</td>
  </tr>
  <tr>
    <td>ComputeAdditionalMetrics_CodonUsage( )</td>
  </tr>
  <tr>
    <td>ConsistencyWithCodonFreq( )</td>
    <td rowspan="2">Internal</td>
  </tr>
  <tr>
    <td>CheckNames_tRNAs( )</td>
  </tr>

  <!-- Codon Freq Per Gene -->
  <tr>
    <td rowspan="7">Codon Freq Per Gene</td>
    <td>CodonFreqPerGene.R</td>
    <td>GetCodonFreq( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td rowspan="6">tTE_HelperFunctions_CodonFreq.R</td>
    <td>ExtractCodons( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>ExtractGenes( )</td>
    <td rowspan="5">Internal</td>
  </tr>
  <tr>
    <td>ExtractSequences( )</td>
  </tr>
  <tr>
    <td>CallingEnsembl( )</td>
  </tr>
  <tr>
    <td>CheckFASTAFormat( )</td>
  </tr>
  <tr>
    <td>FromFASTAtoTable( )</td>
  </tr>

  <!-- Codon Pool Contribution -->
  <tr>
    <td rowspan="4">Pool Contribution</td>
    <td rowspan="4">ExamineCodonPoolContribution.R</td>
    <td>ShowPoolContribution( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>ComputeIndividualGeneCorrelation( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>AnalyzeTopGeneImpact( )</td>
  </tr>
  <tr>
    <td>GetOrCompute( )</td>
  </tr>

  <!-- Amino Acid Usage -->
  <tr>
    <td rowspan="3">Amino Acid Usage</td>
    <td rowspan="3">AAUsage.R</td>
    <td>ComputeAAUsage( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>GroupAA( )</td>
    <td rowspan="2">Internal</td>
  </tr>
  <tr>
    <td>RetrieveAAUsageData( )</td>
  </tr>

  <!-- Theoretical Translation Efficiency -->
  <tr>
    <td rowspan="4">Theoretical Translation Efficiency</td>
    <td rowspan="4">TheoreticalTranslationEfficiency.R</td>
    <td>Compute_tTE( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>ComputeStatisticalSignificance( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>ComputeCorrelation( )</td>
  </tr>
  <tr>
    <td>FilterMatrix( )</td>
  </tr>

  <!-- Abundance analysis -->
  <tr>
    <td rowspan="13">Differential Expression Analysis</td>
    <td rowspan="3">DifferentialExpressionAnalysis.R</td>
    <td>RunDEAnalysis( )</td>
    <td rowspan="3">Visible</td>
  </tr>
  <tr>
    <td>ComputeDEResults( )</td>
  </tr>
  <tr>
    <td>PlotDEResults( )</td>
  </tr>
  <tr>
    <td rowspan="10">tTE_HelperFunctions_DEA.R</td>
    <td>ComputeSizeCorrection( )</td>
    <td rowspan="10">Internal</td>
  </tr>
  <tr>
    <td>ComputeAllPairwiseComp( )</td>
  </tr>
  <tr>
    <td>ComputeDESeq2( )</td>
  </tr>
  <tr>
    <td>TargetedApproach( )</td>
  </tr>
  <tr>
    <td>MakeScatterPlot( )</td>
  </tr>
  <tr>
    <td>RunDimReduct( )</td>
  </tr>
  <tr>
    <td>GenerateVolcanoPlot( )</td>
  </tr>
  <tr>
    <td>HeatmapFontSize( )</td>
  </tr>
  <tr>
    <td>ProduceHeatmapDiffExp( )</td>
  </tr>
  <tr>
    <td>ProduceElbowPlot( )</td>
  </tr>

  <!-- Permutation Test -->
  <tr>
    <td rowspan="2">Permutation Test</td>
    <td rowspan="2">PermutationTest.R</td>
    <td>GetPermutationDist( )</td>
    <td rowspan="2">Visible</td>
  </tr>
  <tr>
    <td>ObtainSignificance( )</td>
  </tr>

  <!-- General Functions -->
  <tr>
    <td rowspan="26">General Functions</td>
    <td rowspan="7">GeneralProcessing.R</td>
    <td>Get_tRNAMatrix( )</td>
    <td rowspan="7">Visible</td>
  </tr>
  <tr>
    <td>Set_tRNACutoff( )</td>
  </tr>
  <tr>
    <td>Filter_tRNACuts( )</td>
  </tr>
  <tr>
    <td>Set_tRNAgenes( )</td>
  </tr>
  <tr>
    <td>TransformFormat( )</td>
  </tr>
  <tr>
    <td>MergeMatrices( )</td>
  </tr>
  <tr>
    <td>GroupConditions( )</td>
  </tr>
  <tr>
    <td rowspan="9">tTE_HelperFunctions_GeneralProcessing.R</td>
    <td>TransformCounts( )</td>
    <td rowspan="9">Internal</td>
  </tr>
  <tr>
    <td>CutoffMatrix( )</td>
  </tr>
  <tr>
    <td>FilteringCutoffs( )</td>
  </tr>
  <tr>
    <td>ComputeCorrelations( )</td>
  </tr>
  <tr>
    <td>Selection_Cutoff( )</td>
  </tr>
  <tr>
    <td>Iterate_tRNACutoff( )</td>
  </tr>
  <tr>
    <td>CorrelationCutoffPlot( )</td>
  </tr>
  <tr>
    <td>SelectionCutoffPlot( )</td>
  </tr>
  <tr>
    <td>Cumulative_SelectionCutoffPlot( )</td>
  </tr>
  <tr>
    <td rowspan="6">tTE_HelperFunctions_GeneralInternal.R</td>
    <td>SelectDefaultData( )</td>
    <td rowspan="6">Internal</td>
  </tr>
  <tr>
    <td>IdentifyInputFormat( )</td>
  </tr>
  <tr>
    <td>FilterByMetadata( )</td>
  </tr>
  <tr>
    <td>CheckCodonFreqTable( )</td>
  </tr>
  <tr>
    <td>CheckDataFrame( )</td>
  </tr>
  <tr>
    <td>CheckGeneAnnotation( )</td>
  </tr>
  <tr>
    <td rowspan="4">tTE_HelperFunctions_GeneralTranslation.R</td>
    <td>FeaturesToAA( )</td>
    <td>Visible</td>
  </tr>
  </tr>
    <td>PerformTranslation( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>RetrieveTranslation( )</td>
  </tr>
  <tr>
    <td>IsEnsemblID( )</td>
  </tr>

  <!-- Plots -->
<tr>
    <td rowspan="22">Plots</td>
    <td rowspan="7">tTEscanPlots.R</td>
    <td>tTE_DistributionPlot( )</td>
    <td rowspan="6">Visible</td>
</tr>
<tr>
    <td>tTE_PlotTargetComparison( )</td>
</tr>
<tr>
    <td>tTE_ProportionPlot( )</td>
</tr>
<tr>
    <td>tTE_ScoresPlot( )</td>
</tr>
<tr>
    <td>tTE_CorrelationPlot( )</td>
</tr>
<tr>
    <td>tTE_PermutationPlot( )</td>
</tr>
<tr>
    <td>CorrelationCutoffPlot( )</td>
    <td>Internal</td>
  </tr>

<tr>
    <td rowspan="15">tTE_HelperFunctions_tTEscanRPlots.R</td>
    <td>SavePlot( )</td>
    <td rowspan="15">Internal</td>
</tr>
<tr>
    <td>GetSafeColorScale( )</td>
</tr>
<tr>
    <td>GetAnnotData( )</td>
</tr>
<tr>
    <td>CheckValueInData( )</td>
</tr>
<tr>
    <td>CheckDataInLongFormat( )</td>
</tr>
<tr>
    <td>GetOutputName( )</td>
</tr>
<tr>
    <td>GenerateDistPlot( )</td>
</tr>
<tr>
    <td>DrawBarCountsPlot( )</td>
</tr>
<tr>
    <td>DrawDonutPlot( )</td>
</tr>
<tr>
    <td>ComputeRings( )</td>
</tr>
<tr>
    <td>DrawRadarPlot( )</td>
</tr>
<tr>
    <td>GenerateProportionPlot( )</td>
</tr>
<tr>
    <td>SignificanceSymbol( )</td>
</tr>
<tr>
    <td>Compute_tTE_Significance( )</td>
</tr>
<tr>
    <td>Compute_Boxplot_Significance( )</td>
</tr>

</table>
