# Contributing to CRISPR Guide RNA Designer

Thank you for your interest in contributing!  
Your help—big or small—improves this tool for everyone.

---

## How to Contribute

**1. Open Issues**
- Found a bug? Have a question or feature request?  
  Please open a [GitHub Issue](../../issues/new) with as much detail as possible.

**2. Submit Pull Requests**
- Fork the repository.
- Create a new branch (`git checkout -b my-feature`).
- Make your changes, add tests if needed.
- Commit and push (`git commit -am "Add cool feature"`).
- Open a Pull Request and describe your changes.

**3. Suggestions & Feedback**
- Feedback on usability, UI, docs, or scientific methods is welcome.
- Feel free to comment on issues or open new ones!

---

## Code Style

- Use clear, readable code and comments.
- Keep functions small and focused.
- Follow [PEP8](https://pep8.org/) style for Python.

---

## App Testing

- If possible, test your changes with `uvicorn main:app --reload` from `crispr_app/` before submitting a PR.
- Make sure the app launches and runs without errors.
- Keep the app deterministic: avoid re-introducing Gemini/OpenAI/LLM explainer dependencies.

---

## Community

- Be respectful and open-minded—this is a friendly project for scientists, students, and the curious.
- Credit will be given for all accepted contributions!

---

## License

By contributing, you agree your code will be released under the [MIT License](LICENSE).

---

Thank you for making CRISPR guide design easier for everyone!  
— Dinesh K
