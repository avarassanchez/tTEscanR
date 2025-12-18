# tTEscanR Function Map <a href="https://your-link-here"><img src="https://github.com/user-attachments/assets/8cca530a-f0a5-4284-bf2e-cf030a8193fa" alt="logo_tTEscanR" align="right" width="120"></a>
<img src="https://img.shields.io/badge/Language-R-blue.svg" style="zoom:100%;" />

<table>
  <tr>
    <th>Section</th>
    <th>File</th>
    <th>File Status</th>
    <th>Function</th>
    <th>Function Status</th>
  </tr>

  <!-- tTEscanRObject -->
  <tr>
    <td rowspan="5">tTEscanRObject</td>
    <td>tTEscanRObject.R</td>
    <td>Main</td>
    <td>Create_tTEscanR_Object( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td rowspan="4">tTE_HelperFunctions_tTEscanRObject.R</td>
    <td rowspan="4">Helper</td>
    <td>Update_tTEscanR_Object( )</td>
    <td>Visible</td>
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
    <td>Main</td>
    <td>Run_tTEscanR_pipeline( )</td>
    <td>Visible</td>
  </tr>

  <!-- Codon Usage -->
  <tr>
    <td rowspan="8">Codon Usage</td>
    <td rowspan="2">CodonUsage.R</td>
    <td rowspan="2">Main</td>
    <td>ComputeCodonUsage( )</td>
    <td rowspan="2">Visible</td>
  </tr>
  <tr>
    <td>ComputeAnticodonUsage( )</td>
  </tr>
  <tr>
    <td rowspan="6">tTE_HelperFunctions_CodonUsage.R</td>
    <td rowspan="6">Helper</td>
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
    <td>Main</td>
    <td>GetCodonFreq( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td rowspan="6">tTE_HelperFunctions_CodonFreq.R</td>
    <td rowspan="6">Helper</td>
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
  <tr>
    <td>ExtractCodons( )</td>
    <td>Visible</td>
  </tr>

  <!-- Codon Pool Contribution -->
  <tr>
    <td rowspan="4">Codon Pool Contribution</td>
    <td>ExamineCodonPoolContribution.R</td>
    <td>Main</td>
    <td>ShowPoolContribution( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td rowspan="3">tTE_HelperFuncitons_PoolContribution.R</td>
    <td rowspan="3">Helper</td>
    <td>ComputeIndividualGeneCorrelation( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>ComputeTopNGenes( )</td>
  </tr>
  <tr>
    <td>ComputeWithoutTopNGenes( )</td>
  </tr>

  <!-- Amino Acid Usage -->
  <tr>
    <td rowspan="3">Amino Acid Usage</td>
    <td>AAUsage.R</td>
    <td>Main</td>
    <td>ComputeAAUsage( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td rowspan="2">tTE_HelperFunctions_AAUsage.R</td>
    <td rowspan="2">Helper</td>
    <td>GroupAA( )</td>
    <td rowspan="2">Internal</td>
  </tr>
  <tr>
    <td>RetrieveAAUsageData( )</td>
  </tr>

  <!-- Theoretical Translation Efficiency -->
  <tr>
    <td rowspan="4">Theoretical Translation Efficiency</td>
    <td>TheoreticalTranslationEfficiency.R</td>
    <td>Main</td>
    <td>Compute_tTE( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td rowspan="3">tTE_HelperFunctions_TranslationEfficiency.R</td>
    <td rowspan="3">Helper</td>
    <td>ComputeStatisticalSignificance( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>ComputeCorrelation( )</td>
  </tr>
  <tr>
    <td>FilterMatrix( )</td>
  </tr>

  <!-- DEA -->
  <tr>
    <td rowspan="13">Differential Expression Analysis</td>
    <td rowspan="4">DifferentialExpressionAnalysis.R</td>
    <td rowspan="4">Main</td>
    <td>ComputeSizeCorrection( )</td>
    <td rowspan="4">Visible</td>
  </tr>
  <tr>
    <td>RunDEAnalysis( )</td>
  </tr>
  <tr>
    <td>ComputeDEResults( )</td>
  </tr>
  <tr>
    <td>PlotDEResults( )</td>
  </tr>
  <tr>
    <td rowspan="9">tTE_HelperFunctions_DEA.R</td>
    <td rowspan="9">Helper</td>
    <td>ComputeDESeq2( )</td>
    <td rowspan="9">Internal</td>
  </tr>
  <tr>
    <td>TargetedApproach( )</td>
  </tr>
  <tr>
    <td>ExploratoryApproach( )</td>
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
    <td rowspan="2">Main</td>
    <td>GetPermutationDist( )</td>
    <td rowspan="2">Visible</td>
  </tr>
  <tr>
    <td>ObtainSignificance( )</td>
  </tr>

  <!-- General Functions -->
  <tr>
    <td rowspan="24">General Functions</td>
    <td rowspan="7">GeneralProcessing.R</td>
    <td rowspan="7">Helper</td>
    <td>Set_tRNACutoff( )</td>
    <td rowspan="7">Visible</td>
  </tr>
  <tr>
    <td>Get_tRNAMatrix( )</td>
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
    <td rowspan="7">tTE_HelperFunctions_GeneralProcessing.R</td>
    <td rowspan="7">Helper</td>
    <td>TransformCounts( )</td>
    <td rowspan="7">Internal</td>
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
    <td>Generate_tRNAMatrix_SE( )</td>
  </tr>
  <tr>
    <td>Generate_tRNAMAtrix( )</td>
  </tr>
  <tr>
    <td rowspan="3">tTE_HelperFunctions_GeneralInternal_InputData.R</td>
    <td rowspan="3">Helper</td>
    <td>SelectDefaultData( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>IdentifyInputFormat( )</td>
  </tr>
  <tr>
    <td>FilterByMetadata( )</td>
  </tr>
  <tr>
    <td rowspan="3">tTE_HelperFunctions_GeneralInternal_Checks.R</td>
    <td rowspan="3">Helper</td>
    <td>CheckCodonFreqTable( )</td>
    <td rowspan="3">Internal</td>
  </tr>
  <tr>
    <td>CheckDataFrame( )</td>
  </tr>
  <tr>
    <td>CheckGeneAnnotation( )</td>
  </tr>
  <tr>
    <td rowspan="4">tTE_HelperFunctions_GeneralTranslation.R</td>
    <td rowspan="4">Helper</td>
    <td>PerformTranslation( )</td>
    <td rowspan="2">Internal</td>
  </tr>
  <tr>
    <td>RetrieveTranslation( )</td>
  </tr>
  <tr>
    <td>FeaturesToAA( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>IsEnsemblID( )</td>
    <td>Internal</td>
  </tr>

  <!-- Plots -->
<tr>
    <td rowspan="19">Plots</td>
    <td rowspan="7">tTEscanPlots.R</td>
    <td rowspan="7">Main</td>
    <td>tTE_DistributionPlot( )</td>
    <td rowspan="6">Visible</td>
</tr>
<tr>
    <td>tTE_CompareTargetToMean( )</td>
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
    <td rowspan="12">tTE_HelperFunctions_tTEscanRPlots.R</td>
    <td rowspan="12">Helper</td>
    <td>SavePlot( )</td>
    <td rowspan="12">Internal</td>
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

</table>
