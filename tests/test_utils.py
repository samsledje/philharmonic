from unittest.mock import MagicMock, patch

import matplotlib.pyplot as plt
import networkx as nx
import pytest

from philharmonic.utils import (
    add_GO_function,
    calculate_graph_triangles,
    clean_top_terms,
    create_rainbow_colorbar,
    download_file_safe,
    filter_proteins_GO,
    get_cluster_top_terms,
    get_node_colors,
    hash_cluster,
    load_cluster_json,
    nx_graph_cluster,
    parse_GO_database,
    parse_GO_graph,
    parse_GO_map,
    plot_cluster,
    plot_degree,
    print_cluster,
    subset_GO_graph,
    write_cluster_cytoscape,
    write_cluster_fasta,
)


# Test print_cluster
def test_print_cluster(capsys):
    cluster = {
        "members": ["protein1", "protein2", "protein3"],
        "GO_terms": {"GO:0001": 3, "GO:0002": 2},
        "graph": [("protein1", "protein2", 0.8), ("protein2", "protein3", 0.7)],
        "recipe": {"degree": {"0.75": ["protein4"]}},
        "llm_name": "Protein Cluster Alpha",
        "llm_explanation": "This cluster contains proteins involved in cellular respiration.",
        "llm_confidence": "Low",
    }
    go_database = {"GO:0001": "Term 1", "GO:0002": "Term 2"}

    def assert_in(s):
        assert "Cluster of 3 proteins" in s
        assert "Top Terms:" in s
        assert "Cluster Name: Protein Cluster Alpha" in s
        assert "Edges: 2" in s
        assert "Triangles: 0" in s
        assert "Max Degree: 2" in s
        assert "This cluster contains proteins involved in cellular respiration." in s

    cluster_string = print_cluster(cluster, go_database, n_terms=2, return_str=True)
    assert_in(cluster_string)

    cluster_string = print_cluster(cluster, go_database, n_terms=2, return_str=False)
    captured = capsys.readouterr()
    assert cluster_string is None
    assert_in(captured.out)


# Test nx_graph_cluster
def test_nx_graph_cluster():
    cluster = {
        "members": ["protein1", "protein2", "protein3"],
        "graph": [("protein1", "protein2", 0.8), ("protein2", "protein3", 0.7)],
        "recipe": {"degree": {"0.75": ["protein4"]}},
    }
    full_G = nx.Graph()
    full_G.add_weighted_edges_from(
        [
            ("protein1", "protein2", 0.8),
            ("protein2", "protein3", 0.7),
            ("protein3", "protein4", 0.6),
        ]
    )

    G = nx_graph_cluster(cluster, full_G=full_G, use_recipe_nodes=True)
    assert len(G.nodes) == 4
    assert len(G.edges) == 3


# Test load_cluster_json
def test_load_cluster_json(tmp_path):
    json_file = tmp_path / "test_clusters.json"
    json_file.write_text('{"cluster1": {"members": ["protein1", "protein2"]}}')
    clusters = load_cluster_json(json_file)
    assert "cluster1" in clusters
    assert clusters["cluster1"]["members"] == ["protein1", "protein2"]


# Test parse_GO_graph
def test_parse_GO_graph(tmp_path):
    go_graph_file = tmp_path / "test_go_graph.obo"
    go_graph_file.write_text(
        "[Term]\n"
        "id: GO:0001\n"
        "name: term1\n"
        "namespace: biological_process\n"
        "is_a: GO:0002\n"
        "\n"
        "[Term]\n"
        "id: GO:0002\n"
        "name: term2\n"
        "namespace: biological_process\n"
    )
    go2children, go2desc = parse_GO_graph(go_graph_file)
    assert go2children == {"GO:0002": ["GO:0001"]}
    assert go2desc == {
        "GO:0001": ("biological_process", "term1"),
        "GO:0002": ("biological_process", "term2"),
    }


# Test subset_GO_graph
def test_subset_GO_graph(tmp_path):
    go_graph_file = tmp_path / "test_go_graph.obo"
    go_graph_file.write_text(
        "[Term]\n"
        "id: GO:0001\n"
        "name: term1\n"
        "namespace: biological_process\n"
        "is_a: GO:0002\n"
        "\n"
        "[Term]\n"
        "id: GO:0002\n"
        "name: term2\n"
        "namespace: biological_process\n"
    )
    result = subset_GO_graph(go_graph_file, ["GO:0002"])
    assert len(result) == 2
    assert ("GO:0002", "biological_process", "term2", "GO:0001") in result
    assert ("GO:0001", "biological_process", "term1", "") in result


# Test parse_GO_database
def test_parse_GO_database(tmp_path):
    go_db_file = tmp_path / "test_go_db.obo"
    go_db_file.write_text(
        "[Term]\nid: GO:0001\nname: term1\n\n[Term]\nid: GO:0002\nname: term2\n"
    )
    go_db = parse_GO_database(go_db_file)
    assert go_db == {"GO:0001": "term1", "GO:0002": "term2"}


# Test parse_GO_map
def test_parse_GO_map(tmp_path):
    go_map_file = tmp_path / "test_go_map.csv"
    go_map_file.write_text(
        "seq,manual_annot,pfam_list,GO_list\n"
        "protein1,,,GO:0007420;GO:0046705\n"
        "protein2,,,GO:0007566;GO:0003256;GO:0042492\n"
    )
    go_map = parse_GO_map(go_map_file)
    assert go_map == {
        "protein1": ["GO:0007420", "GO:0046705"],
        "protein2": ["GO:0007566", "GO:0003256", "GO:0042492"],
    }


# Test plot_degree
def test_plot_degree():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5)])
    plot_degree(G, name="Test Graph")
    assert plt.gcf().number > 0  # Check if a figure was created


# Test plot_cluster
def test_plot_cluster():
    cluster = {
        "members": ["protein1", "protein2", "protein3"],
        "graph": [("protein1", "protein2", 0.8), ("protein2", "protein3", 0.7)],
        "recipe": {"degree": {"0.75": ["protein4"]}},
    }
    full_graph = nx.Graph()
    full_graph.add_weighted_edges_from(
        [
            ("protein1", "protein2", 0.8),
            ("protein2", "protein3", 0.7),
            ("protein3", "protein4", 0.6),
        ]
    )
    plot_cluster(cluster, full_graph, name="Test Cluster", show=False)
    assert plt.gcf().number > 0  # Check if a figure was created


# Test write_cluster_fasta
def test_write_cluster_fasta(tmp_path):
    cluster_file = tmp_path / "test_clusters.json"
    cluster_file.write_text('{"cluster1": {"members": ["protein1", "protein2"]}}')
    sequence_file = tmp_path / "test_sequences.fasta"
    sequence_file.write_text(">protein1\nMKVLWAAS\n>protein2\nMKLPVRGS\n")
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    write_cluster_fasta(cluster_file, sequence_file, directory=output_dir)
    assert (output_dir / "cluster_cluster1.fasta").exists()


# Test write_cluster_cytoscape
def test_write_cluster_cytoscape(tmp_path):
    cluster = {
        "members": ["protein1", "protein2", "protein3"],
        "graph": [("protein1", "protein2", 0.8), ("protein2", "protein3", 0.7)],
        "recipe": {"degree": {"0.75": ["protein4"]}},
    }
    full_G = nx.Graph()
    full_G.add_weighted_edges_from(
        [
            ("protein1", "protein2", 0.8),
            ("protein2", "protein3", 0.7),
            ("protein3", "protein4", 0.6),
        ]
    )
    outfile = tmp_path / "cytoscape_input.txt"
    write_cluster_cytoscape(cluster, full_G, outfile=outfile)
    assert outfile.exists()
    assert (tmp_path / "cytoscape_input.nodes.txt").exists()


# Test create_rainbow_colorbar
def test_create_rainbow_colorbar():
    fig = create_rainbow_colorbar()
    assert isinstance(fig, plt.Figure)


# Test add_GO_function
def test_add_GO_function():
    cluster = {"members": ["protein1", "protein2", "protein3"]}
    go_map = {
        "protein1": ["GO:0001", "GO:0002"],
        "protein2": ["GO:0002", "GO:0003"],
        "protein3": ["GO:0001", "GO:0003"],
    }
    result = add_GO_function(cluster, go_map)
    assert result == {"GO:0001": 2, "GO:0002": 2, "GO:0003": 2}


# Test calculate_graph_triangles
def test_calculate_graph_triangles():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 3)])
    assert calculate_graph_triangles(G) == 2


# Test get_cluster_top_terms
@pytest.mark.parametrize(
    "n_terms,expected",
    [
        (3, [("GO:0001", 5), ("GO:0003", 4), ("GO:0002", 3)]),
        (-1, [("GO:0001", 5), ("GO:0003", 4), ("GO:0002", 3), ("GO:0004", 2)]),
    ],
)
def test_get_cluster_top_terms(n_terms, expected):
    cluster = {"GO_terms": {"GO:0001": 5, "GO:0002": 3, "GO:0003": 4, "GO:0004": 2}}
    result = get_cluster_top_terms(cluster, N=n_terms)
    assert result == expected


# Test hash_cluster
def test_hash_cluster():
    protein_list = ["protein1", "protein2", "protein3"]
    result = hash_cluster(protein_list)
    assert isinstance(result, int)
    assert result == hash_cluster(
        ["protein3", "protein1", "protein2"]
    )  # Order shouldn't matter


# Test filter_proteins_GO
@pytest.mark.parametrize(
    "proteins,expected",
    [
        (["protein1", "protein2", "protein3"], {"protein1", "protein2"}),
        (["protein4", "protein5"], set()),
    ],
)
def test_filter_proteins_GO(proteins, expected, tmp_path):
    go_filter_f = tmp_path / "go_filter.txt"
    go_map_f = tmp_path / "go_map.csv"
    go_database_f = tmp_path / "go_database.obo"

    # Create mock files
    go_filter_f.write_text("GO:0001\nGO:0002")
    go_map_f.write_text(
        "seq,manual_annot,pfam_list,GO_list\n"
        "protein1,,,GO:0001\n"
        "protein2,,,GO:0002\n"
        "protein3,,,GO:0003"
    )
    go_database_f.write_text(
        "[Term]\n"
        "id: GO:0001\n"
        "name: term1\n"
        "namespace: biological_process\n"
        "is_a: GO:0002\n"
        "\n"
        "[Term]\n"
        "id: GO:0002\n"
        "name: term2\n"
        "namespace: biological_process\n"
        "[Term]\n"
        "id: GO:0003\n"
        "name: term3"
        "namespace: biological_process\n"
    )

    result = filter_proteins_GO(proteins, go_filter_f, go_map_f, go_database_f)
    assert result == expected


# Test clean_top_terms
def test_clean_top_terms():
    cluster = {
        "members": ["protein1", "protein2", "protein3"],
        "GO_terms": {"GO:0001": 3, "GO:0002": 2},
    }
    go_db = {"GO:0001": "Term 1", "GO:0002": "Term 2"}
    result = clean_top_terms(cluster, go_db)
    assert result == "Term 1"


# Test get_node_colors
def test_get_node_colors():
    cluster = {
        "members": ["protein1", "protein2"],
        "recipe": {"degree": {"0.75": ["protein3"]}},
    }
    result = get_node_colors(cluster)
    assert result == {"protein1": "blue", "protein2": "blue", "protein3": "red"}


# Test download_file_safe
@patch("requests.get")
def test_download_file_safe(mock_get, tmp_path):
    mock_response = MagicMock()
    mock_response.headers = {"content-length": "12"}
    mock_response.iter_content.return_value = [b"test content"]
    mock_get.return_value = mock_response

    output_path = tmp_path / "test_file.txt"
    result = download_file_safe("http://example.com/file", output_path)

    assert result
    assert output_path.exists()
    assert output_path.read_text() == "test content"
