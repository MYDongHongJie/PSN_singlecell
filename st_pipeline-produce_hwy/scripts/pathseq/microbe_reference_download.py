import requests
import os.path
import os
import re

urls = [
    ('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt', 'bacteria'),
    ('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt', 'fungi'),
    ('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt', 'viral'),
    ('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt', 'archaea')
]

for url, name in urls:
    filename = os.path.join(name, os.path.basename(url))
    response = requests.get(url)
    os.makedirs(name, exist_ok=True)
    with open(filename, 'wb') as f:
        f.write(response.content)

    if name == 'bacteria':
        with open(filename, 'r') as f:
            lines = f.readlines()

        representative_lines = [line for line in lines if 'representative' in line]

        output_filename = os.path.join(name, 'representative.txt')
        with open(output_filename, 'w') as f:
            f.writelines(representative_lines)

        print(f"Saved {len(representative_lines)} representative lines from {filename} to {output_filename}")

        with open(output_filename, 'r') as f:
            lines = f.readlines()

        pattern = r'(https?://\S+)'
        links = []
        for line in lines:
            match = re.findall(pattern, line)
            if match:
                link = match[0] + '/' + match[0].split('/')[-1] + '_genomic.fna.gz'
                links.append(link)

        output_filename = os.path.join(name, 'links.txt')
        with open(output_filename, 'w') as f:
            f.writelines('\n'.join(links))

        print(f"Saved {len(links)} links from {filename} to {output_filename}")
    else:
        with open(filename, 'r') as f:
            lines = f.readlines()

        pattern = r'(https?://\S+)'
        links = []
        for line in lines:
            match = re.findall(pattern, line)
            if match:
                link = match[0] + '/' + match[0].split('/')[-1] + '_genomic.fna.gz'
                links.append(link)

        output_filename = os.path.join(name, 'links.txt')
        with open(output_filename, 'w') as f:
            f.writelines('\n'.join(links))

        print(f"Saved {len(links)} links from {filename} to {output_filename}")

    links_file = os.path.join(name, 'links.txt')
    output_dir = os.path.join(name, 'data')
    os.makedirs(output_dir, exist_ok=True)

    with open(links_file, 'r') as f:
        genome_urls = f.readlines()

    for genome_url in genome_urls:
        genome_url = genome_url.strip()  # remove leading/trailing whitespaces and newline characters
        filename = os.path.basename(genome_url)
        output_path = os.path.join(output_dir, filename)
        response = requests.get(genome_url)
        with open(output_path, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {output_path}")
