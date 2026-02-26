"""
OpenWebUI Tool: Structure Similarity Search

Enables chemists to search for similar molecular structures directly from chat.
Uses the Structure Search API to find molecules with similar ECFP4 fingerprints.

Phase 2, Task 2.4

Usage in chat:
  @tool search_similar_structures
  Query: Find molecules similar to CN1CC[C@]2...
"""

import os
import logging
import requests
from typing import Optional
from pydantic import BaseModel, Field

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Valve(BaseModel):
    """Tool configuration valves (configurable in OpenWebUI admin)."""

    STRUCTURE_SEARCH_API_URL: str = Field(
        default="http://localhost:8000/v1/structure",
        description="Base URL for Structure Search API",
    )
    SIMILARITY_THRESHOLD: float = Field(
        default=0.7,
        ge=0.0,
        le=1.0,
        description="Default similarity threshold (0.0-1.0)",
    )
    MAX_RESULTS: int = Field(
        default=5, ge=1, le=20, description="Maximum number of results to display"
    )
    TIMEOUT_SECONDS: int = Field(
        default=30, ge=5, le=120, description="Request timeout in seconds"
    )


class Tool:
    """
    Structure Similarity Search Tool for OpenWebUI.

    Allows users to search for molecules similar to a query SMILES
    directly from the chat interface.
    """

    val = Valve()

    async def search_similar_structures(
        self,
        query_smiles: str,
        threshold: Optional[float] = None,
        max_results: Optional[int] = None,
    ) -> str:
        """
        Find molecules similar to query structure.

        Args:
            query_smiles: SMILES string of query molecule
            threshold: Optional similarity threshold (overrides default)
            max_results: Optional max results (overrides default)

        Returns:
            Formatted response with search results or error message

        Example:
            >>> await search_similar_structures("CN1CC[C@]2...")
            "🔬 **Similar Structures Found:**

            - **Similarity**: 0.92
              **SMILES**: `CN1CC[C@]2...`
              **Document**: Mesembrine alkaloids...

            - **Similarity**: 0.85
              **SMILES**: `CN1CC[C@@]2...`
              **Document**: Sceletium phytochemistry..."
        """
        api_url = f"{self.val.STRUCTURE_SEARCH_API_URL}/search"
        threshold = threshold or self.val.SIMILARITY_THRESHOLD
        max_results = max_results or self.val.MAX_RESULTS

        try:
            # Call Structure Search API
            response = requests.post(
                api_url,
                json={
                    "query_smiles": query_smiles,
                    "threshold": threshold,
                    "top_k": max_results,
                },
                timeout=self.val.TIMEOUT_SECONDS,
            )

            if response.status_code != 200:
                return f"❌ Structure search failed: {response.text}"

            data = response.json()
            results = data.get("results", [])

            if not results:
                return (
                    f"🔬 **No Similar Structures Found**\n\n"
                    f"No molecules with similarity ≥ {threshold} were found.\n\n"
                    f"**Query SMILES**: `{query_smiles[:60]}...`\n\n"
                    f"Try lowering the threshold or check if the molecule exists in the database."
                )

            # Format results for display
            lines = ["🔬 **Similar Structures Found**:\n"]

            for i, result in enumerate(results, 1):
                similarity = result.get("similarity", 0)
                smiles = result.get("smiles", "")
                doc_title = result.get("document_title", "Unknown")
                metadata = result.get("molecule_metadata", {})

                # Add similarity bar visualization
                bar_length = int(similarity * 20)
                bar = "█" * bar_length + "░" * (20 - bar_length)

                # Build result card
                lines.append(f"**{i}. Similarity: {similarity:.2f}** `{bar}`")
                lines.append(
                    f"   **SMILES**: `{smiles[:60]}{'...' if len(smiles) > 60 else ''}`"
                )
                lines.append(f"   **Document**: {doc_title}")

                # Add metadata if available
                if metadata:
                    props = metadata.get("properties", {})
                    if props:
                        mw = props.get("mw")
                        logp = props.get("logp")
                        if mw:
                            lines.append(f"   **MW**: {mw} Da")
                        if logp:
                            lines.append(f"   **LogP**: {logp}")

                lines.append("")  # Blank line between results

            # Add summary
            lines.append("---")
            lines.append(
                f"**Summary**: Found {len(results)} similar structures "
                f"(threshold ≥ {threshold})"
            )

            return "\n".join(lines)

        except requests.Timeout:
            logger.error(f"Search timeout after {self.val.TIMEOUT_SECONDS}s")
            return (
                f"❌ **Search Timeout**\n\n"
                f"The structure search API did not respond within {self.val.TIMEOUT_SECONDS} seconds. "
                f"Please try again or contact the administrator."
            )
        except requests.ConnectionError as e:
            logger.error(f"Connection error: {e}")
            return (
                f"❌ **Connection Error**\n\n"
                f"Could not connect to Structure Search API at `{api_url}`.\n\n"
                f"Please verify the API is running and the URL is correct in tool settings."
            )
        except Exception as e:
            logger.error(f"Search failed: {e}")
            return f"❌ **Search Failed**: {str(e)}"

    def get_molecule_info(self, doc_id: str) -> str:
        """
        Retrieve molecule information by document ID.

        Args:
            doc_id: Document ID to lookup

        Returns:
            Formatted molecule information or error message
        """
        api_url = f"{self.val.STRUCTURE_SEARCH_API_URL}/{doc_id}"

        try:
            response = requests.get(api_url, timeout=self.val.TIMEOUT_SECONDS)

            if response.status_code == 404:
                return f"❌ **Molecule Not Found**: No molecule data for document `{doc_id}`"

            if response.status_code != 200:
                return f"❌ **Lookup Failed**: {response.text}"

            data = response.json()

            # Format as markdown card
            lines = [
                f"🧪 **Molecule Information**",
                f"",
                f"**Document ID**: `{doc_id}`",
                f"**SMILES**: `{data.get('smiles', 'N/A')}`",
            ]

            if data.get("document_title"):
                lines.append(f"**Document**: {data['document_title']}")

            metadata = data.get("molecule_metadata", {})
            if metadata:
                lines.append("")
                lines.append("**Properties**:")

                props = metadata.get("properties", {})
                if props.get("mw"):
                    lines.append(f"- Molecular Weight: {props['mw']} Da")
                if props.get("logp"):
                    lines.append(f"- LogP: {props['logp']}")
                if props.get("compound_class"):
                    lines.append(f"- Class: {props['compound_class']}")

            return "\n".join(lines)

        except requests.Timeout:
            return f"❌ **Timeout**: Molecule lookup timed out"
        except requests.ConnectionError:
            return f"❌ **Connection Error**: Cannot reach Structure Search API"
        except Exception as e:
            return f"❌ **Error**: {str(e)}"


# =============================================================================
# Tool Registration for OpenWebUI
# =============================================================================


def get_tool_spec():
    """
    Return OpenWebUI tool specification.

    This is automatically discovered by OpenWebUI when the tool is registered.
    """
    return {
        "name": "search_similar_structures",
        "description": "Search for molecules similar to a query SMILES structure",
        "parameters": {
            "type": "object",
            "properties": {
                "query_smiles": {
                    "type": "string",
                    "description": "SMILES string of the query molecule",
                },
                "threshold": {
                    "type": "number",
                    "description": "Similarity threshold (0.0-1.0), default 0.7",
                },
                "max_results": {
                    "type": "integer",
                    "description": "Maximum results to return, default 5",
                },
            },
            "required": ["query_smiles"],
        },
    }


if __name__ == "__main__":
    # Test example
    import asyncio

    tool = Tool()

    # Test with mesembrine
    test_smiles = "CN1CC[C@]2([C@@H]1CC(=O)CC2)C3=CC(=C(C=C3)OC)OC"

    print("Testing Structure Search Tool")
    print("=" * 60)
    print(f"Query SMILES: {test_smiles[:60]}...")
    print("=" * 60)

    result = asyncio.run(tool.search_similar_structures(test_smiles))
    print(result)
