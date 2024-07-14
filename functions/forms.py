from typing import Dict, Union, List



def nullify_empty_string(value:str) -> Union[str, None]:
    """
    This function takes a dictionary and replaces empty strings with None.

    Args:
        value: the value to be checked
    
    Returns:
        str: the value if it is not an empty string
        None: if the value is an empty string
    """
    if value == '':
        return None
    else:
        return value



def get_form_item(request, field:str) -> Union[str, None]:
    """
    This function takes a request object and a field name and returns the value of the field from the form.

    Args:
        request: request object
        field: field name
    Returns:
        str: value of the field
    """
    if field in request.form:
        value = request.form[field]
        return nullify_empty_string(value)
    else:
        return None


def get_querystring_item(request, field:str) -> Union[str, None]:
    """
    This function takes a request object and a field name and returns the value of the field from the querystring.
    """
    if field in request.form:
        value = request.args.get[field]
        return nullify_empty_string(value)
    else:
        return None
    


def get_request_data(request, form:Dict) -> Dict[str, Union[str, None]]:
    """
    This function takes a request object and a list of fields and returns a dictionary of the field values.

    Args:
        request: request object
        fields: list of field names

    Returns:
        dict: dictionary of field values
    """
    data = {}
    if request.method == 'POST':
        field_getter = get_form_item
    elif request.method == 'GET':
        field_getter = get_querystring_item
    else:
        field_getter = None

    if not field_getter:
        return data
    
    for field in form['fields']:
        value = field_getter(request, field)
        if value and 'default_value' in form['fields'][field]:
            if value == form['fields'][field]['default_value']:
                value = None
        data[field] = value
    return data