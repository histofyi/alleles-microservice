from typing import List, Dict, Optional, Callable, Union

from flask import render_template


def render(template_name : str, variables : Dict) -> str:
    """
    This function is used to retrieve the user cookie it does no checking of the privileges of that user

    Args:
        template_name (string): name of the template, with or without the file extension if it's an HTML template
        variables (dictionary): dictionary of variables to be shown in the templated response

    Returns:
        A string containing the rendered html (or other format)
    """
    if ".html" not in template_name:
        template_name += ".html"
    return render_template(template_name, **variables)