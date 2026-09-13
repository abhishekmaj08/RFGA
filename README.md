<h1 align="center"> RFGA: Tumor and Molecular Subtype Mapping in Spatial Transcriptomics</h1>

<p>
<b>RFGA (Random Forest–Genetic Algorithm)</b> is a robust, generalizable supervised machine-learning framework designed for
<b>spot-level tumor identification and molecular subtype mapping in spatial transcriptomics</b>.
By combining <b>Random Forest classification</b> with <b>Genetic Algorithm–based feature optimization</b>,
RFGA provides a flexible framework for generating interpretable spatial predictions of tumor identity and molecular subtype.
</p>

<p>
Importantly, RFGA is designed as a <b>cancer-type-independent framework</b> rather than a tool restricted to a single malignancy.
The current implementation and application of RFGA focus on <b>papillary thyroid carcinoma (PTC)</b>, where it is used to
identify tumor regions and map molecular subtypes at the spatial-spot level.
The framework is intended to be <b>expanded to additional cancer types</b> as appropriate biological reference data,
pathology information, and cancer-specific molecular features are incorporated.
</p>

<h2> Highlights</h2>

<ul>
  <li> <b>Spot-level tumor vs. normal prediction</b></li>
  <li>
     <b>Spatial mapping of PTC thyroid tumor subtypes:</b>
    <ul>
      <li>BRAF-like</li>
      <li>RAS-like</li>
      <li>RET-fusion</li>
      <li>Other subtypes</li>
    </ul>
  </li>
  <li> <b>Tumor probability and boundary mapping</b></li>
  <li> <b>Pathology-guided feature optimization</b></li>
  <li> <b>Interpretable gene signatures</b> for downstream biological analysis</li>
</ul>

<hr>

<h1> Repository Overview</h1>

<p>
The RFGA repository is organized into <b>four main directories</b>:
</p>

<pre><code>RFGA/
│
├── Code/
│   └── R scripts for running RFGA and generating the optimized models
│
├── GA_Optimized_ClassifierModels/
│   └── Genetic Algorithm–optimized tumor classifier models
│
├── AccessoryFiles/
│   └── Files necessary for running the RFGA scripts
│
└── SpatialTranscriptomics_Data_M1_Sample/
    └── Thyroid cancer spatial transcriptomics data for sample M1
</code></pre>

<h3> 1. <code>Code</code></h3>

<p>
Contains the <b>R scripts</b> used to generate final RFGA annotations and to generate the optimized classification models.
</p>

<h3> 2. <code>GA_Optimized_ClassifierModels</code></h3>

<p>
Contains the <b>Genetic Algorithm–optimized tumor classifier models</b> used by RFGA.
</p>

<h3> 3. <code>AccessoryFiles</code></h3>

<p>
Contains the <b>files necessary for running the RFGA scripts</b>.
</p>

<h3> 4. <code>SpatialTranscriptomics_Data_M1_Sample</code></h3>

<p>
Contains our <b>thyroid cancer spatial transcriptomics data for sample M1</b>, which is provided as an example dataset for running RFGA.
</p>

<blockquote>
<b>📥 Large Files:</b> Some files exceed GitHub's 25 MB file-size limit. For these files, a
<b>Dropbox download link</b> is provided in the repository so that the required files can be downloaded separately.
</blockquote>

<hr>

<h1> Generating Final Annotation Results</h1>

<p> <b>RFGA is used to annotate spatial transcriptomics data.</b> In the following sections, we demonstrate the usage of RFGA using our own thyroid cancer spatial transcriptomics data as an example. We also provide instructions for using RFGA to annotate <b>your own spatial transcriptomics data</b>. </p>

<h2>1. Generating annotation for Our Data</h2>

<p>
The script <code>FinalAnnotationPrediction_OurData.R</code> produces the final annotation of a spatial transcriptomics sample using the RFGA model.
We use our thyroid cancer spatial transcriptomics sample <b>M1</b> as the example.
</p>

<h3>Step 1 — Prepare the Required Files</h3>

<p>
Place <b>all required files from the <code>Code</code> folder</b> in the same directory from which the analysis will be run.
</p>

<p>
Use the data in the folder <code>SpatialTranscriptomics_Data_M1_Sample/</code> for our example:
</p>

<h3>Step 2 — Run RFGA</h3>

<p>Open a terminal in the directory containing the required files and run:</p>

<pre><code>Rscript FinalAnnotationPrediction_OurData.R M1 <b>path_to_</b>SpatialTranscriptomics_Data_M1_Sample<b>_directory</b>
</code></pre>

<h3> Expected Output</h3>

<p>The analysis generates <b>three output files</b>:</p>

<h4>1. <code>M1_RFGA_FinalAnnotation.csv</code></h4>

<p>
A CSV file listing the <b>probabilities and final annotation of each spatial spot</b>.
</p>

<h4>2. <code>M1_RFGA_Normal_Tumor_Annotation.tiff</code></h4>

<p>
An image of the spatial sample showing <b>normal and tumor spots</b>.
</p>

<h4>3. <code>M1_RFGA_Tumor_SubType_Annotation.tiff</code></h4>

<p>
An image of the spatial sample showing <b>normal spots and the different tumor-subtype spots</b>.
</p>

<hr>

<h2>2. Generating annotation for Own Data</h2>

<p>
Users can also run RFGA using their <b>own spatial transcriptomics data</b>.
</p>

<h3>Step 1 — Prepare the Required Files</h3>

<p>
Place <b>all required files from the <code>Code</code> folder</b> in the same directory from which the analysis will be run.
</p>

<p>
Instead of the provided M1 data, use your own data. Make sure your spatial data contains the .h5 file.
</p>

<h3>Step 2 — Run RFGA</h3>

<p>Open a terminal in the directory containing the required files and run:</p>

<pre><code>Rscript FinalAnnotationPrediction_OurData.R YourSampleName <b>path_to_</b>Your_Spatial_Sample<b>_directory</b>
</code></pre>

<h3> Expected Output</h3>

<p>The analysis generates <b>three output files</b> same as above, but prefixed by your sample name</p>

<p> For example, if <code>YourSampleName</code> is <code>Sample1</code>, the following files will be generated: </p>
<pre><code>Sample1_RFGA_FinalAnnotation.csv Sample1_RFGA_Normal_Tumor_Annotation.tiff Sample1_RFGA_Tumor_SubType_Annotation.tiff </code></pre>

<p>
The current RFGA implementation and classifier models provided in this repository are based on the
<b>thyroid cancer application</b> described above.
</p>

<hr>

<h1> Running the Genetic Algorithm Scripts</h1>

<p>
The repository also provides the Genetic Algorithm (GA) scripts used to generate the optimized RFGA classifier models.
</p>

<h2>Step 1 — Prepare the Required Files</h2>

<p>
Place <b>all files from the <code>AccessoryFiles</code> folder</b> in the same directory as the GA script that you want to run.
</p>

<h2>Step 2 — Run the Appropriate GA Script</h2>

<p>
Run the required R script from the terminal using:
</p>

<pre><code>Rscript &lt;script_name&gt;.R
</code></pre>

<p>The available GA scripts are:</p>

<h3> <code>NormalTumor_GA.R</code></h3>

<p>
Runs the <b>Genetic Algorithm</b> to generate the optimized model for classifying
<b>Normal and Tumor spatial spots</b>.
</p>

<h3> <code>BrafModel_GA.R</code></h3>

<p>
Runs the <b>Genetic Algorithm</b> to generate the optimized model for classifying
<b>Normal and BRAF-subtype tumor spatial spots</b>.
</p>

<h3> <code>RasModel_GA.R</code></h3>

<p>
Runs the <b>Genetic Algorithm</b> to generate the optimized model for classifying
<b>Normal and RAS-subtype tumor spatial spots</b>.
</p>

<h3> <code>RETModel_GA.R</code></h3>

<p>
Runs the <b>Genetic Algorithm</b> to generate the optimized model for classifying
<b>Normal and RET-subtype tumor spatial spots</b>.
</p>

<hr>
<h2> Cite Our Paper</h2>
<p> If you use <b>RFGA</b>, its code, models, or associated resources in your research, please cite our paper: </p>
<blockquote> Abhishek Majumdar, Yanan Song, Yimin Liu, Matthew D. Ringel, Lang Li, Chongwen Dong, Zhongchao Mai, Wei Xia &amp; Lijun Cheng*, <b>RFGA: Integrating Random Forests and Genetic Algorithms for Tumor and Tumor-Subtype Mapping in Spatial Transcriptomics</b>, <i>under review</i>. </blockquote>
<p> Thank you for supporting our work by citing RFGA in publications and research that use this framework. </p>
