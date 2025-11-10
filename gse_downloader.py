import os
import tarfile
import zipfile
from ftplib import FTP
import pandas as pd

class GSEDownloader:
    """
    Clase para descargar archivos de GEO (GSE) incluyendo:
    - SOFT, MINiML, Series Matrix
    - Archivos suplementarios
    - Matrices crudas de conteos (si existen)
    """

    def __init__(self, gse_id, dest_dir):
        self.gse_id = gse_id
        self.dest_dir = dest_dir
        os.makedirs(dest_dir, exist_ok=True)
        self.downloaded_files = []
        self.raw_count_df = None

    def download(self):
        ftp = FTP("ftp.ncbi.nlm.nih.gov")
        ftp.login()
        base_path = f"/geo/series/{self.gse_id[:-3]}nnn/{self.gse_id}/"

        # -----------------------------
        # 1️⃣ Archivos principales (SOFT, MINiML, Series Matrix)
        # -----------------------------
        ftp.cwd(base_path)
        files = ftp.nlst()
        for filename in files:
            if filename.endswith((".soft.gz", ".xml.gz", "_series_matrix.txt.gz")):
                self._download_file(ftp, filename)

        # -----------------------------
        # 2️⃣ Archivos suplementarios (suppl/)
        # -----------------------------
        try:
            ftp.cwd(base_path + "suppl")
            supp_files = ftp.nlst()
            for filename in supp_files:
                self._download_file(ftp, filename)
        except Exception as e:
            print(f"⚠️ No se encontró carpeta suppl/ o error: {e}")

        ftp.quit()

        # -----------------------------
        # 3️⃣ Buscar matrices crudas
        # -----------------------------
        self._load_raw_counts()
        return self.downloaded_files, self.raw_count_df

    def _download_file(self, ftp, filename):
        local_path = os.path.join(self.dest_dir, filename)
        print(f"📥 Descargando {filename} ...")
        with open(local_path, "wb") as f:
            ftp.retrbinary(f"RETR {filename}", f.write)
        self.downloaded_files.append(local_path)
        print(f"✅ Guardado en {local_path}")

        # Descomprimir automáticamente si es .tar.gz o .zip
        if filename.endswith(".tar.gz"):
            with tarfile.open(local_path, "r:gz") as tar:
                tar.extractall(self.dest_dir)
        elif filename.endswith(".zip"):
            with zipfile.ZipFile(local_path, "r") as zip_ref:
                zip_ref.extractall(self.dest_dir)

    def _load_raw_counts(self):
        for file in os.listdir(self.dest_dir):
            if ("raw" in file.lower() or "count" in file.lower()) and file.endswith((".txt", ".txt.gz", ".csv")):
                raw_path = os.path.join(self.dest_dir, file)
                print(f"📊 Leyendo matriz de conteos desde {file} ...")
                if file.endswith(".gz"):
                    self.raw_count_df = pd.read_csv(raw_path, sep="\t", compression="gzip", index_col=0)
                else:
                    self.raw_count_df = pd.read_csv(raw_path, sep="\t", index_col=0)
                print(f"✅ Matriz de conteos cargada: {self.raw_count_df.shape}")
                break  # Solo la primera matriz cruda encontrada
