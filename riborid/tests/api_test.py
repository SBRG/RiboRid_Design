import os
from riborid import get_auth
import pytest


class TestRequest:
    def __init__(self):
        self.pwd = os.environ['PASSWORD']
        self.usr = 'sapoudel'
        self.token_file = os.environ['CLIENT_INFO']


@pytest.fixture(scope='class')
def create_request(request):
    request.cls.test_api = TestRequest()


@pytest.mark.usefixtures("create_request")
class TestAPI:
    def test_token_file_exists(self):
        assert os.path.isfile(self.test_api.token_file)

    def test_token_file(self):
        with open(self.test_api.token_file, 'r') as f:
            self.client_info = f.readline().strip()
            get_auth.get_authentication(self.test_api.client_info,
                                        self.test_api.usr, self.test_api.pwd)
        assert os.path.isfile('access_token.txt')

    def test_token_content(self):
        with open('access_token.txt', 'r') as f:
            line = f.readline().strip()
            assert line.isalnum()
