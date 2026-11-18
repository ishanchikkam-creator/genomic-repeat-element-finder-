import http.server
import socketserver
import urllib.parse
import repeat_finder as grf

PORT = 8000

# HTML template with basic styling
HTML_PAGE = """
<!DOCTYPE html>
<html>
<head>
    <title>Genomic Repeat Finder</title>
    <style>
        body { font-family: Arial, sans-serif; max-width: 800px; margin: 0 auto; padding: 20px; background: #f4f4f9; }
        .container { background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
        h1 { color: #333; margin-top: 0; }
        .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 15px; margin-bottom: 20px; }
        label { font-weight: bold; font-size: 14px; display: block; margin-bottom: 5px; }
        input[type="number"] { padding: 8px; width: 100%; box-sizing: border-box; border: 1px solid #ccc; border-radius: 4px; }
        textarea { width: 100%; padding: 10px; border: 1px solid #ccc; border-radius: 4px; box-sizing: border-box; font-family: monospace; font-size: 14px; }
        input[type="submit"] { background: #4CAF50; color: white; border: none; padding: 12px 20px; cursor: pointer; border-radius: 4px; font-size: 16px; width: 100%; font-weight: bold; }
        input[type="submit"]:hover { background: #45a049; }
        .results { margin-top: 25px; background: #e8f5e9; padding: 15px; border-radius: 5px; border: 1px solid #c8e6c9; }
        .error { background: #ffebee; border: 1px solid #ffcdd2; color: #c62828; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Genomic Repeat Finder</h1>
        <form method="POST" action="/">
            <div class="grid">
                <div><label>k-mer length:</label><input type="number" name="k_len" value="10"></div>
                <div><label>Min distance:</label><input type="number" name="min_dist" value="50"></div>
                <div><label>Min unit len:</label><input type="number" name="min_unit" value="1"></div>
                <div><label>Max unit len:</label><input type="number" name="max_unit" value="6"></div>
                <div><label>Min repeats:</label><input type="number" name="min_rep" value="2"></div>
            </div>
            
            <label>Paste FASTA text or raw DNA sequence:</label>
            <textarea name="seq" rows="8">>seq1\nATATATATATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCATATATATAT</textarea><br><br>
            
            <input type="submit" value="Run Analysis">
        </form>
        {results}
    </div>
</body>
</html>
"""

class RequestHandler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        # Serve the page when someone visits the URL
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        # This is the line that was fixed (.replace instead of .format)
        html_out = HTML_PAGE.replace("{results}", "")
        self.wfile.write(html_out.encode('utf-8'))

    def do_POST(self):
        # Handle the form submission
        content_length = int(self.headers['Content-Length'])
        post_data = self.rfile.read(content_length).decode('utf-8')
        parsed = urllib.parse.parse_qs(post_data)
        
        try:
            # Extract form values
            k_len = int(parsed.get('k_len', ['10'])[0])
            min_dist = int(parsed.get('min_dist', ['50'])[0])
            min_unit = int(parsed.get('min_unit', ['1'])[0])
            max_unit = int(parsed.get('max_unit', ['6'])[0])
            min_rep = int(parsed.get('min_rep', ['2'])[0])
            seq_input = parsed.get('seq', [''])[0].strip()
            
            if not seq_input:
                raise ValueError("Sequence cannot be empty.")

            # Run the algorithms
            records = grf.parse_sequences(seq_input)
            analyses = []
            for seq_id, seq in records:
                length = len(seq)
                dispersed = grf.find_dispersed_repeats(seq, k_len=k_len, min_dist=min_dist)
                tandem = grf.find_tandem_repeats(seq, min_unit_len=min_unit, max_unit_len=max_unit, min_repeats=min_rep)
                analyses.append({"id": seq_id, "length": length, "dispersed": dispersed, "tandem": tandem})
                
            # Generate output
            gff3_text = grf.to_gff3(analyses, k_len=k_len)
            
            res_html = f'''
            <div class="results">
                <h2>Analysis Complete ✅</h2>
                <label>GFF3 Output:</label>
                <textarea rows="12" readonly>{gff3_text}</textarea>
            </div>
            '''
        except Exception as e:
            res_html = f'<div class="results error"><h2>Error</h2><p>{str(e)}</p></div>'
            
        # Send results back to the browser
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        # This is the line that was fixed (.replace instead of .format)
        html_out = HTML_PAGE.replace("{results}", res_html)
        self.wfile.write(html_out.encode('utf-8'))

if __name__ == "__main__":
    with socketserver.TCPServer(("", PORT), RequestHandler) as httpd:
        print(f"Server running! Open this link in your browser: http://localhost:{PORT}")
        httpd.serve_forever()