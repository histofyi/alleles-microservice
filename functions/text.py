
def slugify(text:str) -> str:
    slug_char = '_'
    to_replace = [' ','-','.',',','[',']','{','}','(',')','/','\\','*',':']
    for replace_me in to_replace:
        if replace_me in text:
            text = text.replace(replace_me, slug_char)
    while slug_char*2 in text:
        text = text.replace(slug_char*2, slug_char)
    if text[0] == slug_char:
        text = text[1:]
    if text[-1] == slug_char:
        text = text[:-1]
    return text.lower()