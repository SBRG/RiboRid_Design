#!/usr/bin/env python

from base64 import b64encode
import json
from urllib import request, parse
import getpass


def get_authentication(client_info, idt_username, idt_password):
    """
    Create the HTTP request and send it, parsing the response for the access token.
    The body_dict will also contain the fields "expires_in" (value of 3600) and "token_type" (value of "Bearer").
    If the request fails for any reason, an exception will be thrown that contains debugging information.
    """
    authorization_string = b64encode(bytes(client_info, "utf-8")).decode()

    request_headers = {"Content-Type": "application/x-www-form-urlencoded",
                       "Authorization": "Basic " + authorization_string}
    data_dict = {"grant_type": "password",
                 "scope": "test",
                 "username": idt_username,
                 "password": idt_password}

    request_data = parse.urlencode(data_dict).encode()
    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token",
                                   data=request_data,
                                   headers=request_headers,
                                   method="POST")

    response = request.urlopen(post_request)
    body = response.read().decode()

    if response.status != 200:
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)

    body_dict = json.loads(body)
    return body_dict["access_token"]


# running from cmd line
if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description='Generate tokens for IDT API. Requires IDT account with '
                                            'API access enabled. For info on how to set up API visit: '
                                            'https://idtdna.com/pages/apidoc')
    p.add_argument('--usr', help='IDT username', required=True)
    p.add_argument('-i', help='Identity file containing a single line with API client_id and client_secret '
                              'separated by \':\'',
                   required=True)

    idt_pwd = getpass.getpass()
    params = vars(p.parse_args())

    with open(params['i'], 'r') as f:
        cinfo = f.readline().strip()

    token = get_authentication(cinfo, params['usr'], idt_pwd)

    with open('access_token.txt', 'w') as out:
        out.write(token)
