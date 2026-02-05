from pathlib import Path


def test_requirements_do_not_include_llm_sdks():
    req = Path('crispr_app/requirements.txt').read_text().lower()
    banned = ['google-generativeai', 'openai', 'anthropic', 'langchain']
    for name in banned:
        assert name not in req
