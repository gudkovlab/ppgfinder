<h2>GitHub Project Description: PPG Finder</h2>

<h3>Overview</h3>
<p>PPG Finder is a toolkit tailored for identifying and analyzing processed pseudogenes in reference genomes. It consists of scripts that aid in genomic data processing and annotation analysis, streamlining the pseudogene identification process.</p>

<h3>Repository Contents</h3>
<ol>
    <li><strong>GTF.py</strong> (Originally authored by Kamil Slowikowski, modified by Gudkov's Lab)
        <ul>
            <li><strong>Description</strong>: Adapted for Python 3, this script is used for parsing and processing GTF files.</li>
            <li><strong>Modifications</strong>: Adapted for Python 3 compatibility.</li>
            <li><strong>License</strong>: Inherits the public domain status from the original script. The original script is at <a href="https://gist.github.com/slowkow/8101481">Kamil Slowikowski's GTF gist</a>.</li>
        </ul>
    </li>
    <li><strong>junction_library_generation.py</strong>
        <ul>
            <li><strong>Description</strong>: This script generates an exon junction library.</li>
            <li><strong>Key Features</strong>: Extraction and processing of exon information from GTF files. Generation of an exon junction library.</li>
            <li><strong>Usage</strong>: 
                <pre>python junction_library_generation.py &lt;gtf_file&gt; &lt;genome_fasta_file&gt; &lt;output_pkl_file&gt; &lt;output_fasta_file&gt; &lt;overhang&gt;</pre>
                <ul>
                    <li>&lt;gtf_file&gt;: Path to the GTF file.</li>
                    <li>&lt;genome_fasta_file&gt;: Path to the genome FASTA file.</li>
                    <li>&lt;output_pkl_file&gt;: Path for the output pickle file.</li>
                    <li>&lt;output_fasta_file&gt;: Path for the output FASTA file.</li>
                    <li>&lt;overhang&gt;: Length parameter for the number of nucleotides in the sequence to the left or right of exact junction.</li>
                </ul>
            </li>
        </ul>
    </li>
    <li><strong>genome_hits_generation.sh</strong>
        <ul>
            <li><strong>Description</strong>: A Bash script for conducting BLAST searches on species genomes, using the exon junction library generated by the <code>junction_library_generation.py</code> script.</li>
            <li><strong>Key Features</strong>: Automated creation and management of BLAST databases. Execution of BLAST searches using the generated exon junction library.</li>
            <li><strong>Usage</strong>:
                <pre>./genome_hits_generation.sh &lt;genome_directory&gt; &lt;junction_library_path&gt; [output_filename]</pre>
                <ul>
                    <li>&lt;genome_directory&gt;: Directory containing the genome files.</li>
                    <li>&lt;junction_library_path&gt;: Path to the exon junction library.</li>
                    <li>[output_filename]: Optional parameter for the name of the output file.</li>
                </ul>
            </li>
        </ul>
    </li>
</ol>

<h3>Installation and Setup</h3>
<p>Clone the repository:</p>
<pre>git clone https://github.com/gudkovlab/ppgfinder.git</pre>
<p><strong>Dependencies</strong>:</p>
<ul>
    <li>Python (3.x+)</li>
    <li>pandas</li>
    <li><a href="https://www.htslib.org/">Samtools</a> (required for <code>junction_library_generation.py</code>)</li>
    <li><a href="https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html">BLAST+ tools</a> (required for <code>genome_hits_generation.sh</code>)</li>
</ul>
<p><strong>Setting Up</strong>:</p>
<ul>
    <li>Install Python and pandas via pip. Install Samtools and BLAST+ from their official sources.</li>
    <li>Make scripts executable using <code>chmod +x script_name</code>.</li>
</ul>

<h3>Usage Workflow</h3>
<ol>
    <li><strong>Create Exon Junction Library</strong>:
        <ul>
            <li>First, use <code>junction_library_generation.py</code> to generate an exon junction library from the GTF and genome FASTA files.</li>
        </ul>
    </li>
    <li><strong>Perform BLAST Search</strong>:
        <ul>
            <li>Then, use <code>genome_hits_generation.sh</code> to perform BLAST searches on the species genome using the exon junction library.</li>
        </ul>
    </li>
</ol>

<h3>Contributing</h3>
<p>Contributions to PPG Finder can be made by forking the repository and submitting a pull request.</p>

<h3>License</h3>
<p>PPG Finder, excluding <code>GTF.py</code>, is released under the MIT License. <code>GTF.py</code> is a modified version, inheriting the public domain status from the original script. The original script is available at <a href="https://gist.github.com/slowkow/8101481">Kamil Slowikowski's GTF gist</a>.</p>

<h3>Contact</h3>
<p>For queries or feedback, please open an issue in this GitHub repository.</p>

<h3>Acknowledgments</h3>
<p>We thank Kamil Slowikowski for the original <code>GTF.py</code> script.</p>
