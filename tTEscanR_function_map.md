## tTEscanR Function Map

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
    <td>DefineNewData( )</td>
    <td>Internal</td>
  </tr>
  <tr>
    <td>Update_tTEscanR_Object( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>CheckInputCombinations( )</td>
    <td>Internal</td>
  </tr>
  <tr>
    <td>IsIn_tTEscanR_Object( )</td>
    <td>Internal</td>
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
    <td rowspan="7">Codon Usage</td>
    <td rowspan="2">CodonUsage.R</td>
    <td rowspan="2">Main</td>
    <td>ComputeCodonUsage( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>ComputeAnticodonUsage( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td rowspan="5">tTE_HelperFunctions_CodonUsage.R</td>
    <td rowspan="5">Helper</td>
    <td>ComputeCodonExonicBackground( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>ComputeCorrelationBackground( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>ComputeMeanUsage( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>ComputeAdditionalMetrics_CodonUsage( )</td>
    <td>Visible</td>
  </tr>
  <tr>
    <td>ConsistencyWithCodonFreq( )</td>
    <td>Internal</td>
  </tr>

  <!-- Codon Freq Per Gene -->
  <tr>
    <td rowspan="8">Codon Freq Per Gene</td>
    <td>Co
