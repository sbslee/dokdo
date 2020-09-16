import nbformat as nbf

def create_report(**kwargs):
    nb = nbf.v4.new_notebook()

    m_intro = """\
This is an auto-generated notebook."""

    c_import = """\
from qiime2 import Visualization"""

    c_table_qzv = """\
Visualization.load('table.qzv')"""

    c_rep_seqs_qzv = """\
Visualization.load('rep-seqs.qzv')"""

    c_alpha_rarefaction_qzv = """\
Visualization.load('alpha-rarefaction.qzv')"""

    c_shannon_group_significance_qzv = """\
Visualization.load('core-metrics-results/shannon_group-significance.qzv')"""

    c_bray_curtis_emperor_qzv = """\
Visualization.load('core-metrics-results/bray_curtis_emperor.qzv')"""

    c_taxa_bar_plots_qzv = """\
Visualization.load('taxa-bar-plots.qzv')"""

    code = """\
%pylab inline
hist(normal(size=2000), bins=50);"""

    nb['cells'] = [nbf.v4.new_markdown_cell(m_intro),
                   nbf.v4.new_code_cell(c_import),
                   nbf.v4.new_code_cell(c_table_qzv),
                   nbf.v4.new_code_cell(c_rep_seqs_qzv),
                   nbf.v4.new_code_cell(c_alpha_rarefaction_qzv),
                   nbf.v4.new_code_cell(c_shannon_group_significance_qzv),
                   nbf.v4.new_code_cell(c_bray_curtis_emperor_qzv),
                   nbf.v4.new_code_cell(c_taxa_bar_plots_qzv),
                   nbf.v4.new_code_cell(code)]

    nbf.write(nb, kwargs['o_path'])
