""" Functions for generating the html report at the end of the iBench analysis pipeline.
"""
import os
import webbrowser

from ibench.constants import ENDC_TEXT, OKCYAN_TEXT

def create_html_report(config, figures):
    """ Function to create the final html report and open it in the brower.

    Parameters
    ----------
    config : inspire.config.Config
        The Config object for the whole pipeline.
    figures : dict
        A dictionary containing all of the required plots.
    """
    html_string = ('''
    <html>
        <head>
            <link 
                rel="stylesheet"
                href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css"
            >
            <style>
                body{
                    margin:0 100;
                    padding-bottom: 50px;
                    background:whitesmoke;
                }
                h2{
                    color: firebrick;
                    text-align: center;
                }
                h3{
                    color: firebrick;
                    text-align: center;
                }
                h4{
                    color: firebrick;
                    text-align: center;
                }
                h2:after
                {
                    content:' ';
                    display: block;
                    border:0.25px solid #696969;
                    position: absolute;
                    width: 60%;
                    margin-top: 2px;
                    left: 20%;
                }
            </style>
        </head>
        <body>
            <h2>
                iBench Performance Report for ''' + config.identifier + '''</h2>
            <h3>
                Precision-Recall Curves
            </h3>
            <p>
                This plots the precision (correct identifications divided by all identifications)
                against the recall (true identifications divided by possible identifications)
                for each identification method as we vary the scoring the threshold from the
                maximum to the minimum score in the data.
            </p>
            <center>
    ''' + figures['pr'] +
    '''
            </center>
            <h3>
                Receiver-Operator Curves
            </h3>
            <p>
                This plots the true positive rate (correct identifications divided by correct
                identifications and missed identifications) against the false positive rate
                (incorrect identifications divided by incorrect identifications and unidentifiable
                scan for which no peptide was identified) for each identification method as
                we vary the scoring the threshold from the maximum to the minimum score
                in the data.
            </p>
            <center>
    ''' + figures['roc'] +
    '''
            </center>
            <h3>
                Score Distribution
            </h3>
            <p>
                This plots distributions of scores for correct, incorrect, and decoy PSMs for
                each stratum.
            </p>
            <center>
    ''' + figures['distro']
    )
    if 'fdr' in figures:
        html_string += (
            '''
                </center>
                <h3>
                    False Discovery Rate Estimation
                </h3>
                <p>
                    This plots the estimated FDR against observed FDR for each method.
                </p>
                <center>
            ''' + figures['fdr']
        )

    html_string += ('''
            </center>
            <h3>
                Confounding Variables
            </h3>
            <p>
                This plots search engine score distributions against the values of possible
                confounding variables for correct and incorrect PSMs for each identification
                method.
            </p>
    ''')


    for name, fig in figures['conf'].items():
        html_string += ('''
                    <h4>
            ''' + name +
            '''
                    </h4>
                    <center>
            ''' + fig + '</center>'
        )


    html_string += ('''
        </body>
    </html>
    ''')

    output_path = f'{config.output_folder}/ibench-report.html'
    with open(output_path, 'w', encoding='UTF-8') as output_file:
        output_file.write(html_string)

    print(
        OKCYAN_TEXT +
        '\tReport generated.' +
        ENDC_TEXT
    )

    webbrowser.open(
        'file://' + os.path.realpath(output_path),
        new=2
    )
