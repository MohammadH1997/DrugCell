import json

import requests

AUTHORIZATION_TOKEN = ''
COSMIC_SESSION = ''
SESSION_ID = ''
PAGE_SMITH = ''


def x(gene_name):
    url = f"https://cancer.sanger.ac.uk/cmc/backend/api/visualisation/gene/{gene_name}/"

    payload = {}
    headers = {
        'accept': 'application/json, text/plain, */*',
        'accept-language': 'en-US,en;q=0.9,fa;q=0.8',
        'authorization': f'Bearer {AUTHORIZATION_TOKEN}',
        # 'cookie': '',
        'user-agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0.0.0 Safari/537.36'
    }

    response = requests.request(method="GET",
                                url=url,
                                headers=headers,
                                data=payload,
                                cookies={
                                    'Pagesmith': PAGE_SMITH,
                                    'cosmic_session': COSMIC_SESSION,
                                    # 'CookieControl': '{"necessaryCookies":["cosmic_session","genome_version","ordering-*","visibility-*","_hjIncludedInSample","Pagesmith","piwik_ignore"],"optionalCookies":{},"initialState":{"type":"closed"},"statement":{},"consentDate":1731584416661,"consentExpiry":180,"interactedWith":true,"user":"XXXX"}',
                                    'sessionid': SESSION_ID,
                                })
    # print(response.text)
    response.raise_for_status()

    data = response.json()
    references = data.get('bokeh', dict()).get('doc', dict()).get('roots', dict()).get('references', list())
    cds = [reference for reference in references if reference["type"] == "ColumnDataSource"][0]
    cds_data = cds.get('attributes', dict()).get('data', dict())
    # print(cds_data)
    aa_pos = cds_data.get('aa_pos', list())
    # shield_col_bar_fill_colour = cds_data.get('shield_col_bar_fill_colour', list())
    shield_col_raw = cds_data.get('shield_col_raw', list())

    result = dict()
    for pos, tier in zip(aa_pos, shield_col_raw):
        if tier.lower() == 'other':
            continue
        n = int(tier.lower().replace('tier ', '')) if tier.lower() != 'other' else 0
        result[pos] = n
    return result


def main():
    gene_list = []
    with open("data/gene_API.txt", "r") as file:
        for line in file:
            gene_list.append(line.strip())
    try:
        with open("data/codon_list.txt", "r") as file:
            result = json.load(file)
    except Exception as ex:
        result = dict()
    for gene in gene_list:
        if gene in result:
            print(f'skipping {gene}')
            continue
        print(gene)
        try:
            result[gene] = x(gene)
            with open("data/codon_list.txt", "w") as file:
                json.dump(result, file)
        except Exception as ex:
            print(ex)


if __name__ == '__main__':
    main()
