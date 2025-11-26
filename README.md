<h1>BulkDecon</h1>

<p><strong>BulkDecon</strong> is an R package for performing <strong>bulk RNA-seq deconvolution</strong> using <strong>custom single-cell reference profiles</strong>. 
It extends and adapts the core logic from the SpatialDecon framework for more general bulk deconvolution workflows.</p>

<hr>

<h2>🚧 Current Status</h2>

<p>BulkDecon is under active development. Progress so far:</p>

<ul>
  <li>✔️ Package structure scaffolded</li>
  <li>✔️ Core deconvolution logic implemented (LNR model & outlier filtering)</li>
  <li>✔️ Key functionality ported and adapted from SpatialDecon</li>
  <li>✔️ First working commit</li>
</ul>

<p>If you are interested in testing or contributing, please contact <strong>XX</strong>.</p>

<hr>

<h2>📦 Installation</h2>

<pre><code># install.packages("devtools")
devtools::install_github("shaniAmare/BulkDecon", force = TRUE)
</code></pre>

<hr>

<h2>📘 Overview</h2>

<p>BulkDecon provides:</p>

<ul>
  <li>Gene-aligned single-cell reference integration</li>
  <li>Log-normal regression (LNR) deconvolution</li>
  <li>Outlier detection and gene filtering</li>
  <li>Optional TIL-residual weighting</li>
  <li>Custom plotting utilities (e.g., TIL barplots, floret diagrams)</li>
</ul>

<p>A complete vignette and workflow examples will be added once function development is finalized.</p>

<hr>

<h2>📄 License</h2>

<p>MIT</p>

<hr>

<h2>🤝 Contributing</h2>

<p>Contributions are welcome.  
Please open an issue or contact <strong>shani.amarasinghe@monash.edu</strong> if you would like to test or extend BulkDecon.</p>
