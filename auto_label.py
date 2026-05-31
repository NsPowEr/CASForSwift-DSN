import sys, json, os, re
from pathlib import Path
from graphify.build import build_from_json
from graphify.cluster import score_all
from graphify.analyze import god_nodes, surprising_connections, suggest_questions
from graphify.report import generate

def get_label(nodes):
    if not nodes: return "Empty"
    # Try to find common prefix
    parts_list = [re.split(r'[_.:/]', n) for n in nodes[:10]] # look at top 10
    common = parts_list[0]
    for parts in parts_list[1:]:
        new_common = []
        for a, b in zip(common, parts):
            if a.lower() == b.lower(): new_common.append(a)
            else: break
        common = new_common
    
    if common and len(" ".join(common)) > 3:
        label = " ".join(common).title()
        return label
    
    # Fallback: take the most common first part
    first_parts = [p[0] for p in parts_list if p]
    if first_parts:
        most_common = max(set(first_parts), key=first_parts.count)
        return most_common.title()
        
    return nodes[0][:30]

graph_data = json.loads(Path('graphify-out/graph.json').read_text(encoding="utf-8"))
detection  = json.loads(Path('graphify-out/.graphify_detect.json').read_text(encoding="utf-8"))
analysis   = json.loads(Path('graphify-out/.graphify_analysis.json').read_text(encoding="utf-8"))

G = build_from_json(graph_data)
communities = {int(k): v for k, v in analysis['communities'].items()}
cohesion = {int(k): v for k, v in analysis['cohesion'].items()}
tokens = {'input': 0, 'output': 0}

# Generate labels via heuristic
labels = {cid: get_label(nodes) for cid, nodes in communities.items()}

# Regenerate questions
questions = suggest_questions(G, communities, labels)

# Regenerate report
report = generate(G, communities, cohesion, labels, analysis['gods'], analysis['surprises'], detection, tokens, '.', suggested_questions=questions)

Path('graphify-out/GRAPH_REPORT.md').write_text(report, encoding="utf-8")
Path('graphify-out/.graphify_labels.json').write_text(json.dumps({str(k): v for k, v in labels.items()}, ensure_ascii=False), encoding="utf-8")

print("Labels generated and report updated successfully.")
