import json

from selenium import webdriver
from selenium.webdriver.common.by import By
import time
from pathlib import Path

from bs4 import BeautifulSoup

current_dir = Path(__file__).parent.absolute()


USERNAME = "mohammad.h.beigi@ut.ac.ir"
PASSWORD = "8810383770M@"

def get_driver():
    driver = webdriver.Chrome()
    return driver


def login(driver):
    url = f"https://cancer.sanger.ac.uk/cmc/gene/BRAF/codon/603"
    driver.get(url)
    time.sleep(2)
    driver.find_element(by=By.ID, value="id_username").send_keys(USERNAME)
    driver.find_element(by=By.ID, value="id_password").send_keys(PASSWORD)
    driver.find_element(by=By.CSS_SELECTOR, value="input.submit").click()


# gene_name = "BRAF"
# codon_number = 600
def y(driver, gene_name, codon_number, index):
    # from redis import Redis
    # key = str(gene_name) + ":" + str(codon_number)
    # redis_conn = Redis(host="localhost", port=6379)
    # if redis_conn.get(key):
    #     return
    if (gene_name, codon_number) not in {("BIRC3", "533"), ("BIRC3", "536"), ("BIRC3", "491"), ("BIRC3", "525"), ("BIRC3", "532"), ("BIRC3", "514"), ("BIRC3", "494"), ("BIRC3", "541"), ("BIRC3", "526"), ("BIRC3", "545"), ("BIRC3", "472"), ("BIRC3", "537"), ("BIRC3", "513"), ("BIRC3", "540"), ("BIRC3", "498")}:
        return
    url = f"https://cancer.sanger.ac.uk/cmc/gene/{gene_name}/codon/{codon_number}"
    driver.get(url)
    time.sleep(8)

    html_source = driver.page_source
    soup = BeautifulSoup(html_source, "html.parser")

    # gene_dir = current_dir / gene_name
    # gene_dir.mkdir(exist_ok=True)
    # with open(gene_dir / f"{codon_number}.html", "w", encoding="utf-8") as file:
    #     file.write(html_source)

    parent = soup.find("div", class_="grid-generator-parent")
    boxes = parent.find_all("div", class_="grid-overlay") if parent else []
    box_results = []
    for box in boxes:
        chart_name = dict(enumerate(map(lambda x: x.text, box.find("span", class_="chart-name").find_all("h6"))))
        chart_type = box.find("h6", recursive=False).text

        box_style = box["style"]
        if "rgb(143, 25, 27)" in box_style:
            tier = 1
        elif "rgb(211, 87, 39)" in box_style:
            tier = 2
        elif "rgb(218, 217, 130)" in box_style:
            tier = 3
        elif "rgb(253, 251, 250)" in box_style:
            tier = 4
        else:
            tier = None
        _data = {
            "box_type": chart_type,
            "box_name": chart_name,
            "tier": tier
        }
        box_results.append(_data)
    try:
        with open(current_dir / f"data/cosmic_tier_{index}.json", "r") as f:
            result = json.load(f)
    except Exception as e:
        result = dict()
    if gene_name not in result:
        result[gene_name] = dict()
    if codon_number not in result[gene_name]:
        result[gene_name][codon_number] = dict()
    result[gene_name][codon_number] = box_results
    with open(current_dir / f"data/cosmic_tier_{index}.json", "w") as f:
        json.dump(result, f)
    # redis_conn.set(key, 1, ex=60 * 86400)


def chunk_main(data, index=0):
    driver = get_driver()
    login(driver=driver)

    for gene_name, codon_number in data:
        y(driver, gene_name, codon_number, index)
    driver.quit()


def main():
    from multiprocessing import Process
    with open("/home/mohammad/Public/Drug-cell/DrugCell-public/codon_list_filtered.json", "r") as file:
        data = json.load(file)
    chunk_size = 3420
    pairs = list()
    for gene_name, body in data.items():
        for codon_number, _ in body.items():
            pairs.append((gene_name, codon_number))
    chunks = [pairs[i * chunk_size:(i + 1) * chunk_size] for i in range(0, len(pairs) // chunk_size + 1)]
    from pprint import pprint
    pprint(len(chunks))
    processes = []
    for index, chunk in enumerate(chunks):
        processes.append(Process(target=chunk_main, args=(chunk, index,), daemon=True))
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    # chunk_main(data)


if __name__ == "__main__":
    main()
