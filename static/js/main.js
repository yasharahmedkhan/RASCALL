document.addEventListener('DOMContentLoaded', (event) => {
    const searchForm = document.getElementById('search-form');
    const plotForm = document.getElementById('plot-form');

    if (searchForm) {
        searchForm.addEventListener('submit', function(e) {
            e.preventDefault();
            const query = document.getElementById('search-query').value;
            fetch(`/search_results?query=${encodeURIComponent(query)}`)
                .then(response => response.json())
                .then(data => {
                    const resultsContainer = document.getElementById('search-results');
                    resultsContainer.innerHTML = '';
                    data.forEach(result => {
                        const resultElement = document.createElement('div');
                        resultElement.innerHTML = `
                            <h3>${result['IUPAC Name']}</h3>
                            <p>SMILES: ${result['SMILES']}</p>
                            <p>Molecular Formula: ${result['Molecular Formula']}</p>
                            <p>Molecular Weight: ${result['Molecular Weight']}</p>
                        `;
                        resultsContainer.appendChild(resultElement);
                    });
                })
                .catch(error => console.error('Error:', error));
        });
    }

    if (plotForm) {
        plotForm.addEventListener('submit', function(e) {
            e.preventDefault();
            const mol = document.getElementById('mol').value;
            const fg = document.getElementById('fg').value;
            const mf = document.getElementById('mf').value;

            fetch('/plot_rascall', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({ mol, fg, mf }),
            })
            .then(response => response.json())
            .then(item => {
                Bokeh.embed.embed_item(item, "plot-container");
            })
            .catch((error) => {
                console.error('Error:', error);
            });
        });
    }
});