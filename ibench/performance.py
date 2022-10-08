""" Functions for performance analysis of an identification method or methods on iBench
    ground truth datasets.
"""
from math import ceil, floor

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ibench.constants import CONFOUNDING_FEATURE_NAMES, RESIDUE_WEIGHTS, TRANSPLICED_KEY
from ibench.html_report import create_html_report

def _generate_cut_offs(query_table, score_key, n_steps):
    """ Helper function to generate many score cut offs to for creation of precision recall curve.

    Parameters
    ----------
    query_table : pd.DataFrame
        A DataFrame of PSMs.
    score_key : str
        The name of the column of scores for which we are generating a precision-recall curve.

    Returns
    -------
    score_cuts : list of float
        A list of scoring thresholds.
    """
    car_score_min = query_table[score_key].min()
    car_score_max = query_table[score_key].max()
    car_step_size = (car_score_max - car_score_min)/n_steps
    return [car_score_min + (car_step_size*idx) for idx in range(n_steps)]

def create_pr_curve(qt_df, result_dict, stratum):
    """ Function to create a pr curve for a given result.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    result_dict : dictionary
        Dictionary of meta-data about the method being plotted.
    stratum : str
        The stratum for which the ROC curve is being plotted.

    Returns
    -------
    pr_curve : go.Scatter
        Plotly graph object for the ROC curve.
    min_p : float
        The minimum precison value observed for any threshold.
    """
    name = result_dict['name']
    cut_offs = _generate_cut_offs(qt_df, f'{name}Score', 1_000)
    true_count = qt_df[qt_df['trueStratum'] == stratum].shape[0]
    precisions = []
    recalls = []

    for cut_off in cut_offs:
        pred_count, correct_count, _ = _get_counts(
            qt_df, cut_off, name, stratum
        )
        if pred_count > 0:
            precisions.append(correct_count/pred_count)
            recalls.append(correct_count/true_count)
    min_p = min(precisions)
    return go.Scatter(
        x=recalls,
        y=precisions,
        name=f'{name} PR',
        line={'color': result_dict['colour']}
    ), min_p

def create_roc_curve(qt_df, result_dict, stratum):
    """ Function to create a roc curve for a given result.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    result_dict : dictionary
        Dictionary of meta-data about the method being plotted.
    stratum : str
        The stratum for which the ROC curve is being plotted.

    Returns
    -------
    roc_curve : go.Scatter
        Plotly graph object for the ROC curve.
    """
    name = result_dict['name']
    cut_offs = _generate_cut_offs(qt_df, f'{name}Score', 50)
    true_count = qt_df[qt_df['trueStratum'] == stratum].shape[0]
    tprs = []
    fprs = []
    for cut_off in cut_offs:
        pred_count, correct_count, true_negatives = _get_counts(
            qt_df, cut_off, name, stratum, with_true_negatives=True
        )
        false_positives = pred_count - correct_count
        total_negative = true_negatives + false_positives
        if true_count > 0 and total_negative > 0:
            tprs.append(correct_count/true_count)
            fprs.append(false_positives/total_negative)
    max_tpr = max(tprs)
    return go.Scatter(
        x=fprs,
        y=tprs,
        name=f'{name} ROC',
        line={'color': result_dict['colour']}
    ), max_tpr


def plot_precision_recall(qt_df, config):
    """ Function to plot the precision recall curve of the identification method.


    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    config : ibench.config.Config
        The Config object used to run the experiment.
    """
    strata = [x for x in qt_df['trueStratum'].unique() if x != TRANSPLICED_KEY]

    pr_curves = {}
    min_pr = 1.0
    for stratum in strata:
        pr_curves[stratum] = []
        for results in config.benchmark_results:
            pr_trace, pr_val = create_pr_curve(qt_df, results, stratum)
            pr_curves[stratum].append(pr_trace)
            min_pr = min(min_pr, pr_val)
    y_lower_lim = floor(min_pr*10)/10
    fig = make_subplots(
        shared_yaxes=True,
        rows=1,
        cols=len(strata),
        subplot_titles=[x.title() for x in strata],
    )

    for idx, stratum in enumerate(strata):
        for pr_trace in pr_curves[stratum]:
            fig.add_trace(
                pr_trace,
                row=1,
                col=idx+1,
            )

    fig.update_layout(
        width=880*len(strata),
        height=500,
        title_x=0.5,
        plot_bgcolor='rgba(0,0,0,0)',
        yaxis_showticklabels=True,
    )
    for idx in range(2, len(strata)+1):
        fig.update_layout({f'yaxis{idx}_showticklabels': True})
    fig.update_xaxes(
        showline=True,
        linewidth=0.5,
        linecolor='black',
        showgrid=False,
        ticks="outside",
        range=[0,1],
    )
    fig.update_yaxes(
        showline=True,
        linewidth=0.5,
        linecolor='black',
        showgrid=False,
        ticks="outside",
        range=[y_lower_lim, 1.001]
    )
    for i in range(1, len(strata)+1):
        fig['layout'][f'xaxis{i}']['title'] = 'Recall'
        fig['layout'][f'yaxis{i}']['title'] = 'Precision'
    return fig.to_html()

def create_fdr_curve(qt_df, results, stratum):
    """ Function to create a true vs estimated fdr curve for a given result.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    result_dict : dictionary
        Dictionary of meta-data about the method being plotted.
    stratum : str
        The stratum for which the curve is being plotted.

    Returns
    -------
    roc_curve : go.Scatter
        Plotly graph object for the true vs. estimated curve.
    """
    fdrs = []
    name = results['name']
    goal_fdrs = [0.5*i for i in range(1, 11)]
    found_fdrs = []
    for gofdr, fdr_cut in zip(goal_fdrs, [0.005*i for i in range(1, 11)]):
        filtered_df = qt_df[
            (qt_df[f'{name}qValue'] < fdr_cut) &
            (qt_df[f'{name}Stratum'] == stratum)
        ]
        n_found = filtered_df.shape[0]
        correct_count = filtered_df[filtered_df.apply(
            lambda x : x['truePeptide'].replace('I', 'L') == x[f'{name}Peptide'].replace('I', 'L'),
            axis=1
        )].shape[0]
        if n_found > 0:
            fdrs.append(100*(1-(correct_count/n_found)))
            found_fdrs.append(gofdr)
    return go.Scatter(
        x=found_fdrs,
        y=fdrs,
        name=f'{name} FDR',
        line={'color': results['colour']}
    )

def plot_fdr_estimation(qt_df, config):
    """ Function to plot the accuracy of
    """
    strata = [x for x in qt_df['trueStratum'].unique() if x != TRANSPLICED_KEY]
    fdr_curves = {}

    for stratum in strata:
        fdr_curves[stratum] = [go.Scatter(
            x=[i*0.5 for i in range(11)],
            y=[i*0.5 for i in range(11)],
            name='Perfect FDR Estimation',
            line={'color': 'black', 'dash': 'dash'}
        )]
        for results in config.benchmark_results:
            if results['searchEngine'] in ('mascot', 'percolator'):
                fdr_trace = create_fdr_curve(qt_df, results, stratum)
                fdr_curves[stratum].append(fdr_trace)

    fig = make_subplots(
        shared_yaxes=True,
        rows=1,
        cols=len(strata),
        subplot_titles=[x.title() for x in strata],
    )


    for idx, stratum in enumerate(strata):
        for fdr_trace in fdr_curves[stratum]:
            fig.add_trace(
                fdr_trace,
                row=1,
                col=idx+1,
            )

    fig.update_layout(
        width=880*len(strata),
        height=500,
        title_x=0.5,
        plot_bgcolor='rgba(0,0,0,0)',
        yaxis_showticklabels=True,
    )
    for idx in range(2, len(strata)+1):
        fig.update_layout({f'yaxis{idx}_showticklabels': True})
    fig.update_xaxes(
        showline=True,
        linewidth=0.5,
        linecolor='black',
        showgrid=False,
        ticks="outside",
    )
    fig.update_yaxes(
        showline=True,
        linewidth=0.5,
        linecolor='black',
        showgrid=False,
        ticks="outside",
    )
    for i in range(1, len(strata)+1):
        fig['layout'][f'xaxis{i}']['title'] = 'Estimated FDR'
        fig['layout'][f'yaxis{i}']['title'] = 'Observed FDR'

    return fig.to_html()

def analyse_performance(config):
    """ Function to analyse performance of an identification method or methods on
        iBench ground truth datasets.

    Parameters
    ----------
    config : ibench.config.Config
        The Config object used to run the experiment.
    """
    qt_df = pd.read_csv(f'{config.output_folder}/queryTable.csv')
    figures = {}
    figures['pr'] = plot_precision_recall(qt_df, config)
    figures['roc'] = plot_roc(qt_df, config)
    figures['conf'] = plot_confounding(qt_df, config)
    figures['fdr'] = plot_fdr_estimation(qt_df, config)
    figures['distro'] = plot_overall_distro(qt_df, config)
    create_html_report(config, figures)

def plot_cts(qt_df, result_dict, stratum, feature):
    """ Function to create a plot distributions of scores against continuous variables.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    result_dict : dictionary
        Dictionary of meta-data about the method being plotted.
    stratum : str
        The stratum for which the curve is being plotted.
    feature : str
        The feature being plotted.

    Returns
    -------
    traces : list of go.Scatter
        list of Plotly graph objects for the scatter plots.
    """
    name = result_dict['name']
    qt_df = qt_df[qt_df['trueStratum'] == stratum]
    correct_df = qt_df[
        qt_df.apply(
            lambda x : (
                isinstance(x[f'{name}Peptide'], str) and
                x[f'{name}Peptide'].replace('I', 'L') == x['truePeptide'].replace('I', 'L')
            ),
            axis=1
        )
    ]
    incorrect_df = qt_df[
        qt_df.apply(
            lambda x : (
                isinstance(x[f'{name}Peptide'], str) and
                x[f'{name}Peptide'].replace('I', 'L') != x['truePeptide'].replace('I', 'L')
            ),
            axis=1
        )
    ]
    traces = [
        go.Scatter(
            x=correct_df[feature],
            y=correct_df[f'{name}Score'],
            marker={'color': 'ForestGreen'},
            mode='markers',
        ),
        go.Scatter(
            x=incorrect_df[feature],
            y=incorrect_df[f'{name}Score'],
            marker={'color': 'OrangeRed'},
            mode='markers',
        ),
    ]
    if f'{name}DecoyPeptide' in qt_df.columns:
        decoy_df = qt_df[
            qt_df.apply(
                lambda x : (
                    isinstance(x[f'{name}DecoyPeptide'], str)
                ),
                axis=1
            )
        ]
        traces.append(go.Scatter(
            x=decoy_df[feature],
            y=decoy_df[f'{name}DecoyScore'],
            marker={'color': 'DarkMagenta'},
            mode='markers',
        ))

    return traces

def plot_overall_distro(qt_df, config):
    """ Function to plot engine scores for correct and and incorrect PSMs against confounding
        variables.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    config : ibench.config.Config
        The Config object used to run the experiment.
    """
    strata = [x for x in qt_df['trueStratum'].unique() if x != TRANSPLICED_KEY]

    fig = make_subplots(
        rows=1,
        cols=len(config.benchmark_results),
        subplot_titles=[method['name'] for method in config.benchmark_results],
    )
    for idx, results in enumerate(config.benchmark_results):
        for trace in plot_score_distro(qt_df, results):
            if idx > 0:
                trace.update(showlegend=False)
            fig.add_trace(
                trace,
                row=1,
                col=1+idx,
            )
    fig.update_layout(violinmode='group')
    fig.update_layout(
        title='Score Distributions Per Method',
        width=1000,
        height=400*len(strata),
        title_x=0.5,
        plot_bgcolor='rgba(0,0,0,0)',
        yaxis_showticklabels=True,
        # showlegend=False,
    )
    fig.update_xaxes(
        showline=True,
        linewidth=0.5,
        linecolor='black',
        showgrid=False,
        ticks="outside",
        # range=[0,1],
    )

    fig.update_yaxes(
        showline=True,
        linewidth=0.5,
        linecolor='black',
        showgrid=False,
        ticks="outside",
        # range=[0, y_upper_lim]
    )
    fig['layout']['yaxis1']['title'] = 'PSM Score'
    for i in range(1, len(config.benchmark_results)+1):
        fig['layout'][f'xaxis{i}']['title'] = 'Stratum'
    return fig.to_html()

def plot_score_distro(qt_df, results):
    """ Function to create a plot distributions of scores against discrete variables.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    result_dict : dictionary
        Dictionary of meta-data about the method being plotted.
    stratum : str
        The stratum for which the curve is being plotted.

    Returns
    -------
    traces : list of go.Scatter
        list of Plotly graph objects for the violin plots.
    """
    name = results['name']

    qt_df['correct'] = qt_df.apply(
        lambda x : (
            'correct' if isinstance(x[f'{name}Peptide'], str) and
            x[f'{name}Peptide'].replace('I', 'L') == x['truePeptide'].replace('I', 'L')
            else 'incorrect'
        ),
        axis=1
    )
    if f'{name}DecoyScore' in qt_df.columns:
        qt_df = qt_df.apply(
            lambda x : _add_decoys(x, name), axis=1
        )

    fig = px.violin(
        x=qt_df[f'{name}Stratum'],
        y=qt_df[f'{name}Score'],
        color=qt_df['correct'],
        color_discrete_map={
            'correct': 'ForestGreen',
            'incorrect': 'OrangeRed',
            'decoy': 'DarkMagenta'
        },
    )

    return fig['data']


def plot_discrete(qt_df, results, stratum, feature):
    """ Function to create a plot distributions of scores against discrete variables.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    result_dict : dictionary
        Dictionary of meta-data about the method being plotted.
    stratum : str
        The stratum for which the curve is being plotted.
    feature : str
        The feature being plotted.

    Returns
    -------
    traces : list of go.Scatter
        list of Plotly graph objects for the violin plots.
    """
    name = results['name']
    grouped_df = qt_df.groupby(feature, as_index=False)[f'{name}Score'].count()
    grouped_df.columns = [feature, 'count']
    grouped_df = grouped_df[grouped_df['count'] > qt_df.shape[0]/100]
    plot_vals = grouped_df[feature].tolist()

    qt_df = qt_df[qt_df['trueStratum'] == stratum]
    qt_df = qt_df[qt_df[feature].apply(lambda x : x in plot_vals)]

    qt_df['correct'] = qt_df.apply(
        lambda x : (
            'correct' if isinstance(x[f'{name}Peptide'], str) and
            x[f'{name}Peptide'].replace('I', 'L') == x['truePeptide'].replace('I', 'L')
            else 'incorrect'
        ),
        axis=1
    )
    if f'{name}DecoyScore' in qt_df.columns:
        qt_df = qt_df.apply(
            lambda x : _add_decoys(x, name), axis=1
        )

    fig = px.violin(
        x=qt_df[feature],
        y=qt_df[f'{name}Score'],
        color=qt_df['correct'],
        color_discrete_map={
            'correct': 'ForestGreen',
            'incorrect': 'OrangeRed',
            'decoy': 'DarkMagenta',
        },
    )

    return fig['data']

def _add_decoys(df_row, name):
    """ Function to add decoy scores to the main target scores for plotting distributions.

    Parameters
    ----------
    df_row : pd.Series
        A query table row to be modified.
    name : str
        The name of the method.

    Returns
    -------
    df_row : pd.Series
        The updated DataFrame row.
    """
    if np.isnan(df_row[f'{name}Score']) and df_row[f'{name}DecoyScore'] is not None:
        df_row[f'{name}Stratum'] = df_row[f'{name}DecoyStratum']
        df_row[f'{name}Score'] = df_row[f'{name}DecoyScore']
        df_row['correct'] = 'decoy'

    return df_row

def plot_confounding(qt_df, config):
    """ Function to plot engine scores for correct and and incorrect PSMs against confounding
        variables.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    config : ibench.config.Config
        The Config object used to run the experiment.
    """
    strata = [x for x in qt_df['trueStratum'].unique() if x != TRANSPLICED_KEY]

    qt_df['hydrophobicity'] = qt_df['truePeptide'].apply(
        lambda x : ProteinAnalysis(x).gravy()
    )
    qt_df['peptideLength'] = qt_df['truePeptide'].apply(len)
    qt_df['peptideMass'] = qt_df['truePeptide'].apply(
        lambda pep : sum((RESIDUE_WEIGHTS[res] for res in pep))
    )
    all_figs = {}
    for results in config.benchmark_results:
        fig = make_subplots(
            rows=len(strata)*3,
            cols=2,
            row_titles=strata*3,
        )
        for idx, stratum in enumerate(strata):
            pep_mass_traces = plot_cts(qt_df, results, stratum, 'peptideMass')
            hydrophobicity_traces = plot_cts(qt_df, results, stratum, 'hydrophobicity')
            pep_length_traces = plot_discrete(qt_df, results, stratum, 'peptideLength')
            charge_traces = plot_discrete(qt_df, results, stratum, 'charge')
            coverage_traces = plot_cts(qt_df, results, stratum, 'ms2Coverage')
            sig2no_traces = plot_cts(qt_df, results, stratum, 'signalToNoise')
            for trace in pep_length_traces:
                if idx > 0:
                    trace.update(showlegend=False)
                fig.add_trace(
                    trace,
                    row=1+(idx*3),
                    col=1,
                )
            for trace in charge_traces:
                trace.update(showlegend=False)
                fig.add_trace(
                    trace,
                    row=1+(idx*3),
                    col=2,
                )
            for trace in pep_mass_traces:
                trace.update(showlegend=False)
                fig.add_trace(
                    trace,
                    row=2+(idx*3),
                    col=1,
                )
            for trace in hydrophobicity_traces:
                trace.update(showlegend=False)
                fig.add_trace(
                    trace,
                    row=2+(idx*3),
                    col=2,
                )
            for trace in coverage_traces:
                trace.update(showlegend=False)
                fig.add_trace(
                    trace,
                    row=3+(idx*3),
                    col=1,
                )
            for trace in sig2no_traces:
                trace.update(showlegend=False)
                fig.add_trace(
                    trace,
                    row=3+(idx*3),
                    col=2,
                )
        fig.update_layout(
            width=1000,
            height=800*len(strata),
            title_x=0.5,
            plot_bgcolor='rgba(0,0,0,0)',
            yaxis_showticklabels=True,
        )
        fig.update_xaxes(
            showline=True,
            linewidth=0.5,
            linecolor='black',
            showgrid=False,
            ticks="outside",
            # range=[0,1],
        )

        fig.update_yaxes(
            showline=True,
            linewidth=0.5,
            linecolor='black',
            showgrid=False,
            ticks="outside",
            # range=[0, y_upper_lim]
        )
        fig.update_layout(violinmode='group')
        for i in range(len(strata)*len(CONFOUNDING_FEATURE_NAMES)):
            fig['layout'][f'xaxis{i+1}']['title'] = CONFOUNDING_FEATURE_NAMES[
                i%len(CONFOUNDING_FEATURE_NAMES)
            ]
            if i%2 == 0:
                fig['layout'][f'yaxis{i+1}']['title'] = 'Search Engine Score'

        all_figs[results['name']] = fig.to_html()
    return all_figs

def _get_counts(all_df, score_cut_off, name, acc_grp, with_true_negatives=False):
    """ Function to get the predicted and correct counts above a scoring threshold.

    Parameters
    ----------
    all_df : pd.DataFrame
        A DataFrame of all psms.
    score_cut_off : float
        A cut off on psms percolator score.
    score_key : str
        The column name of the percolator scores.
    acc_grp : str
        The accession group we are interested in.

    Returns
    -------
    pred_count : int
        The number of PSMs for the accession group above the threshold.
    correct_count : int
        The number of correct PSMs for the accession group above the threshold.
    """
    filtered_df = all_df[
        (all_df[f'{name}Score'] >= score_cut_off) &
        (all_df[f'{name}Stratum'] == acc_grp)
    ]
    pred_count = filtered_df.shape[0]

    correct_count = filtered_df[filtered_df.apply(
        lambda x : x['truePeptide'].replace('I', 'L') == x[f'{name}Peptide'].replace('I', 'L'),
        axis=1
    )].shape[0]

    if with_true_negatives:
        tn_df = all_df[
            (
                (all_df[f'{name}Score'] < score_cut_off) |
                (all_df[f'{name}Score'].isna())
            )
            & (all_df['trueStratum'] == TRANSPLICED_KEY)
        ]
        return pred_count, correct_count, tn_df.shape[0]

    return pred_count, correct_count, None

def plot_roc(qt_df, config):
    """ Function to plot the receiver operator curve of the identification method.

    Parameters
    ----------
    qt_df : pd.DataFrame
        The Query Table containing ground truth assignments joined to the user's assignments.
    config : ibench.config.Config
        The Config object used to run the experiment.
    """
    strata = [x for x in qt_df['trueStratum'].unique() if x != TRANSPLICED_KEY]
    roc_curves = {}
    max_tpr = 0.0
    for stratum in strata:
        roc_curves[stratum] = []
        for results in config.benchmark_results:
            roc_trace, tpr = create_roc_curve(qt_df, results, stratum)
            roc_curves[stratum].append(roc_trace)
            if tpr > max_tpr:
                max_tpr = tpr

    y_upper_lim = ceil(max_tpr*10)/10
    fig = make_subplots(
        shared_yaxes=True,
        rows=1,
        cols=len(strata),
        subplot_titles=[x.title() for x in strata],
    )

    for idx, stratum in enumerate(strata):
        for pr_trace in roc_curves[stratum]:
            fig.add_trace(
                pr_trace,
                row=1,
                col=idx+1,
            )

    fig.update_layout(
        width=880*len(strata),
        height=500,
        title_x=0.5,
        plot_bgcolor='rgba(0,0,0,0)',
        yaxis_showticklabels=True,
    )
    for idx in range(2, len(strata)+1):
        fig.update_layout({f'yaxis{idx}_showticklabels': True})

    fig.update_xaxes(
        showline=True,
        linewidth=0.5,
        linecolor='black',
        showgrid=False,
        ticks="outside",
        # range=[0,1],
    )

    fig.update_yaxes(
        showline=True,
        linewidth=0.5,
        linecolor='black',
        showgrid=False,
        ticks="outside",
        range=[0, y_upper_lim]
    )
    for i in range(1, len(strata)+1):
        fig['layout'][f'xaxis{i}']['title'] = 'False Positive Rate'
        fig['layout'][f'yaxis{i}']['title'] = 'True Positive Rate'

    return fig.to_html()
