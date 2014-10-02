<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="/">
#ifndef _COMMON_XML_H_
#define _COMMON_XML_H_
<xsl:for-each select="lbe3d">
<xsl:for-each select="variable">
extern <xsl:value-of select="@type"/><xsl:text disable-output-escaping="yes"> </xsl:text><xsl:value-of select="@name"/>;
</xsl:for-each>
</xsl:for-each>
<xsl:for-each select="lbe3d/core">
#ifdef <xsl:value-of select="@name"/>
<xsl:text disable-output-escaping="yes">
</xsl:text>
<xsl:for-each select="variable">
extern <xsl:value-of select="@type"/><xsl:text disable-output-escaping="yes"> </xsl:text><xsl:value-of select="@name"/>;
</xsl:for-each>
#endif</xsl:for-each>
#endif</xsl:template>
</xsl:stylesheet>
