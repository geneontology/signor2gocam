import requests
import os
from requests.models import Response


class SignorDownloader:
    @staticmethod
    def write_response_content(response: Response, file_basename, destination_dir=None):
        if destination_dir is None:
            destination_dir = "resources"
        file_target = os.path.join(destination_dir, file_basename)
        with open(file_target, "wb") as ft:
            ft.write(response.content)
        return file_target

    @staticmethod
    def download_complexes(destination_dir=None):
        response = requests.post(url="https://signor.uniroma2.it/download_complexes.php",
                                 data={"submit": "Download complex data"})
        file_basename = "SIGNOR_complexes.csv"
        return SignorDownloader.write_response_content(response, file_basename, destination_dir=destination_dir)

    @staticmethod
    def download_families(destination_dir=None):
        response = requests.post(url="https://signor.uniroma2.it/download_complexes.php",
                                 data={"submit": "Download protein family data"})
        file_basename = "SIGNOR_PF.csv"
        return SignorDownloader.write_response_content(response, file_basename, destination_dir=destination_dir)
