import os
from riborid import get_auth
import pytest
import requests

class TestRequest:
    def __init__(self):
        self.pwd = os.environ['PASSWORD']
        self.usr = os.environ['USR']
        self.client_info = os.environ['CLIENT_INFO']


@pytest.fixture(scope='class')
def create_request(request):
    request.cls.test_api = TestRequest()
    yield
    os.remove('access_token.txt')

@pytest.mark.usefixtures("create_request")
class TestAPI:
    def test_token_file_exists(self):
        assert os.path.isfile(self.test_api.client_info)

    def test_token_file(self):
        with open(self.test_api.client_info, 'r') as f:
            cinfo = f.readline().strip()
            get_auth.get_authentication(cinfo,
                                        self.test_api.usr, self.test_api.pwd)
        assert os.path.isfile('access_token.txt')

    def test_token_content(self):
        with open('access_token.txt', 'r') as f:
            line = f.readline().strip()
            assert line.isalnum()
